cap mata mata drop OSQP_workspace_abridged()
cap mata mata drop OSQP_csc_matrix()
cap mata mata drop OSQP_setup()
cap mata mata drop OSQP_solve()
cap mata mata drop OSQP_cleanup()
cap mata mata drop OSQP_csc_matrix()
cap mata mata drop OSQP_vec_export()
cap mata mata drop OSQP_int_export()
cap mata mata drop OSQP_csc_export()
cap mata mata drop OSQP_csc_convert()

mata
struct OSQP_workspace_abridged {
    string scalar info_status
    real scalar info_obj_val
    real vector solution_x
    real scalar rc
}

struct OSQP_csc_matrix {
    real vector data
    real vector indices
    real vector indptr
}

void OSQP_setup (string scalar fname,
                 struct OSQP_csc_matrix P,
                 real vector q,
                 struct OSQP_csc_matrix A,
                 real vector u,
                 real vector l)
{
    real scalar fh
    fh = fopen(fname, "rw")
    OSQP_csc_export(fh, P)
    OSQP_csc_export(fh, A)
    OSQP_vec_export(fh, q)
    OSQP_vec_export(fh, u)
    OSQP_vec_export(fh, l)
    fclose(fh)
}

struct OSQP_workspace_abridged scalar OSQP_solve (string scalar fname) {
    struct OSQP_workspace_abridged scalar work
    stata(sprintf(`"plugin call honestosqp_plugin, `"%s"'"', fname))
    OSQP_cleanup(work, fname)
    return(work)
}

void OSQP_cleanup (struct OSQP_workspace_abridged scalar work, string scalar fname) {
    real scalar fh
    colvector C
    fh = fopen(fname, "r")
    C = bufio()
    work.rc = fbufget(C, fh, "%4bu", 1)
    if ( work.rc == 0 ) {
        work.info_status = strtrim(fbufget(C, fh, "%32s", 1))
        if ( work.info_status == "solved" ) {
            work.info_obj_val = fbufget(C, fh, "%8z", 1)
            work.solution_x   = fbufget(C, fh, "%8z", fbufget(C, fh, "%4bu", 1))
        }
    }
    fclose(fh)
    unlink(fname)
}

struct OSQP_csc_matrix scalar OSQP_csc_convert(real matrix A) {
    struct OSQP_csc_matrix scalar S
    real scalar ix, i, j, nnz

    ix  = 0
    nnz = sum(A :!= 0)
    S.data    = J(1, nnz, .)
    S.indices = J(1, nnz, .)
    S.indptr  = J(1, cols(A) + 1, 0)
    for (j = 1; j <= cols(A); j++) {
        for (i = 1; i <= rows(A); i++) {
            if (A[i, j] == 0) continue
            S.data[++ix]  = A[i, j]
            S.indices[ix] = i - 1
        }
        S.indptr[j + 1] = ix
    }

    return(S)
}

void function OSQP_int_export(real scalar fh, real scalar n)
{
    colvector C
    C = bufio()
    fbufput(C, fh, "%4bu", n)
}

void function OSQP_vec_export(real scalar fh, real vector v)
{
    colvector C
    C = bufio()
    fbufput(C, fh, "%4bu", length(v))
    fbufput(C, fh, "%8z",  v)
}

void function OSQP_csc_export(real scalar fh, struct OSQP_csc_matrix scalar A)
{
    colvector C
    C = bufio()
    fbufput(C, fh, "%4bu", length(A.data))
    fbufput(C, fh, "%4bu", length(A.indices))
    fbufput(C, fh, "%4bu", length(A.indptr))
    fbufput(C, fh, "%8z",  A.data)
    fbufput(C, fh, "%4bu", A.indices)
    fbufput(C, fh, "%4bu", A.indptr)
}

end
