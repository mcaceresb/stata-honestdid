cap mata mata drop OSQP()
cap mata mata drop OSQP_get()
cap mata mata drop OSQP_setup()
cap mata mata drop OSQP_solve()
cap mata mata drop OSQP_read()
cap mata mata drop OSQP_cleanup()
cap mata mata drop OSQP_csc_convert()
cap mata mata drop OSQP_vec_export()
cap mata mata drop OSQP_int_export()
cap mata mata drop OSQP_csc_export()
cap mata mata drop OSQP_csc_matrix()
cap mata mata drop OSQP_workspace_abridged()

* Operator Splitting Quadratic Program
* ------------------------------------
*
* Mata interface to OSQP. This solves problems of the form
*
* min_x
*     0.5 x' P x + q' x
* s.t.
*     l <= Ax <= u
*
* with
*
* - P a symmetric positive semi-definite matrix
* - q, x conformable vectors
* - A any matrix conformable with x
* - l, u bounds conformable with Ax
*
* Following this notation, with P, q, A, u, l the corresponding stata
* matrices, run
*
*     result = OSQP(P, q, A, u, l, 1)
*
* result will be a structure of type OSQP_workspace_abridged with components
*
* - rc: Return code for the function (0 means success)
* - info_status: Status of the optimization ('solved' means a solution was found)
* - info_obj_val: Objective value
* - solution_x: Minimizer
*
* Notes
* -----
*
* - For equality constraints, set the corresponding l and u entries to
*   the same number
* - For one-sided inequality constraints, set the corresponding entry
*   of l or u to missing (.); this is interpreted as -infty or + infty,
*   respectively.
*
* Sparse matrices
* ---------------
*
* If you read the OSQP documentation you will note sparse matrices are
* required.  OSQP converts stata matrices to Compressed Sparse Column
* format internally using the OSQP_csc_matrix structure. If you need
* to create your own sparce matrices using the OSQP_csc_matrix structure,
* you can access the internals as follows:
*
*     bufvar = OSQP_setup(P, q, A, u, l)
*     result = OSQP_solve(bufvar, 1)
*
* Where P, A are already of type OSQP_csc_matrix. IMPORTANT: In
* this case, P needs to be upper-triangular. Since P is meant to be
* symmetric, only one triangle is required, and OSQP's internals assume
* you are only passing the upper triangle of the matrix.
*
* Reference
* ---------
*
* See https://osqp.org/docs

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

real matrix OSQP_get(struct OSQP_workspace_abridged scalar work, string scalar element) {
    if ( element == "info_obj_val" ) return(work.info_obj_val);
    if ( element == "solution_x"   ) return(work.solution_x);
    if ( element == "rc"           ) return(work.rc);
}

struct OSQP_workspace_abridged scalar OSQP(real matrix P,
                                           real vector q,
                                           real matrix A,
                                           real vector u,
                                           real vector l,
                                           | real scalar cleanup,
                                           real scalar verbose)
{
    if ( args() < 6 ) cleanup = 0
    if ( args() < 7 ) verbose = 0
    string scalar bufvar
    bufvar = OSQP_setup(OSQP_csc_convert(uppertriangle(P)), q, OSQP_csc_convert(A), u, l, verbose)
    return(OSQP_solve(bufvar, cleanup))
}

string scalar OSQP_setup (struct OSQP_csc_matrix P,
                          real vector q,
                          struct OSQP_csc_matrix A,
                          real vector u,
                          real vector l,
                          | real scalar verbose)
{
    if ( args() < 6 ) verbose = 0

    string scalar buf, bufvar
    real scalar off

    st_numscalar("__honestosqp_rc", .)
    st_numscalar("__honestosqp_obj", .)
    st_numscalar("__honestosqp_max_iter", length(st_numscalar("__honestosqp_max_iter"))? st_numscalar("__honestosqp_max_iter"): .)
    st_matrix("__honestosqp_x", J(1, length(q), .))

    off = 0
    buf = 32 * char(0)

    OSQP_int_export(buf, verbose, off)
    OSQP_csc_export(buf, P, off)
    OSQP_csc_export(buf, A, off)
    OSQP_vec_export(buf, q, off)
    OSQP_vec_export(buf, u, off)
    OSQP_vec_export(buf, l, off)

    if ( st_nobs() < 1 ) stata("qui set obs 1")
    if ( st_global("__HONESTBUFVAR") == "" ) {
        bufvar = st_tempname()
        (void) st_addvar("strL", bufvar)
        st_global("__HONESTBUFVAR", bufvar)
    }
    else {
        bufvar = st_global("__HONESTBUFVAR")
    }
    (void) st_sstore(1, bufvar, buf)
    return(bufvar)
}

struct OSQP_workspace_abridged scalar OSQP_solve(string scalar bufvar, | real scalar cleanup) {
    if ( args() < 2 ) cleanup = 0
    struct OSQP_workspace_abridged scalar work
    stata(sprintf(`"plugin call honestosqp_plugin %s, _plugin_run"', bufvar))
    OSQP_read(work, cleanup)
    return(work)
}

void OSQP_read(struct OSQP_workspace_abridged scalar work, | real scalar cleanup) {
    if ( args() < 2 ) cleanup = 0
    work.rc = st_numscalar("__honestosqp_rc")
    if ( work.rc == 0 ) {
        work.info_status = st_local("__honestosqp_status")
        if ( work.info_status == "solved" ) {
            work.info_obj_val = st_numscalar("__honestosqp_obj")
            work.solution_x   = st_matrix("__honestosqp_x")
        }
    }
    if ( cleanup ) OSQP_cleanup()
}

void OSQP_cleanup() {
    real scalar ix
    ix = _st_varindex(st_global("__HONESTBUFVAR"))
    if ( !missing(ix) ) (void) st_dropvar(ix)
    st_global("__HONESTBUFVAR", "")
    st_local("__honestosqp_status", "")
    st_numscalar("__honestosqp_rc",       J(0, 0, .))
    st_numscalar("__honestosqp_obj",      J(0, 0, .))
    st_numscalar("__honestosqp_max_iter", J(0, 0, .))
    st_matrix("__honestosqp_x", J(0, 0, .))
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

void function OSQP_int_export(string scalar buf, real scalar n, real scalar off)
{
    colvector C
    C = bufio()
    buf = buf + 4 * char(0);
    bufput(C, buf, off, "%4bu", n); off = off + 4;
}

void function OSQP_vec_export(string scalar buf, real vector v, real scalar off)
{
    colvector C
    C = bufio()
    buf = buf + (4 + 8 * length(v)) * char(0);
    bufput(C, buf, off, "%4bu", length(v)); off = off + 4;
    bufput(C, buf, off, "%8z",  v);         off = off + 8 * length(v);
}

void function OSQP_csc_export(string scalar buf, struct OSQP_csc_matrix scalar A, real scalar off)
{
    colvector C
    C = bufio()
    buf = buf + (12 + 8 * length(A.data) + 4 * length(A.indices) + 4 * length(A.indptr)) * char(0);
    bufput(C, buf, off, "%4bu", length(A.data));    off = off + 4;
    bufput(C, buf, off, "%4bu", length(A.indices)); off = off + 4;
    bufput(C, buf, off, "%4bu", length(A.indptr));  off = off + 4;
    bufput(C, buf, off, "%8z",  A.data);            off = off + 8 * length(A.data);
    bufput(C, buf, off, "%4bu", A.indices);         off = off + 4 * length(A.indices);
    bufput(C, buf, off, "%4bu", A.indptr);          off = off + 4 * length(A.indptr);
}
end
