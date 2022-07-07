cap mata mata drop ECOS_workspace_abridged()
cap mata mata drop ECOS_setup()
cap mata mata drop ECOS_solve()
cap mata mata drop ECOS_cleanup()
cap mata mata drop ECOS_code()

* n      = Number of variables
* m      = Number of constraints (from G)
* p      = Number of equality constraints (from A)
* l      = Dimension of positive ortant (for 1st order cone in G)
* ncones = Number of second-order cones (second order cones in G)
* q      = Array of length ncones; q[i] is dimension of ith 2nd-order cone in G
* e      = Number of exponential cones
* G      = csc matrix G
* A      = csc matrix A
* c      = vector for objective (length n)
* h      = vector for conic constraints (length m)
* b      = vector for equality constraint (length p)

mata
struct ECOS_workspace_abridged {
    string scalar info_status
    real scalar info_exitcode
    real scalar info_obj_val
    real vector solution_x
    real scalar rc
}

void ECOS_setup (string scalar fname,
                 real scalar n,
                 real scalar m,
                 real scalar p,
                 real scalar l,
                 real scalar ncones,
                 real vector q,
                 real scalar e,
                 struct OSQP_csc_matrix G,
                 struct OSQP_csc_matrix A,
                 real vector c,
                 real vector h,
                 real vector b)
{
    real scalar fh
    fh = fopen(fname, "rw")
    OSQP_int_export(fh, n)
    OSQP_int_export(fh, m)
    OSQP_int_export(fh, p)
    OSQP_int_export(fh, l)
    OSQP_int_export(fh, ncones)
    ECOS_vec_export(fh, q)
    OSQP_int_export(fh, e)
    OSQP_csc_export(fh, G)
    OSQP_csc_export(fh, A)
    OSQP_vec_export(fh, c)
    OSQP_vec_export(fh, h)
    OSQP_vec_export(fh, b)
    fclose(fh)
}

struct ECOS_workspace_abridged scalar ECOS_solve (string scalar fname) {
    struct ECOS_workspace_abridged scalar work
    stata(sprintf(`"plugin call honestecos_plugin, `"%s"'"', fname))
    ECOS_cleanup(work, fname)
    return(work)
}

void ECOS_cleanup (struct ECOS_workspace_abridged scalar work, string scalar fname) {
    real scalar fh
    colvector C
    fh = fopen(fname, "r")
    C = bufio()
    work.rc = fbufget(C, fh, "%4bu", 1)
    if ( work.rc == 0 ) {
        work.info_exitcode = fbufget(C, fh, "%4bu", 1)
        work.info_status   = ECOS_code(work.info_exitcode)
        work.info_obj_val  = fbufget(C, fh, "%8z", 1)
        work.solution_x    = fbufget(C, fh, "%8z", fbufget(C, fh, "%4bu", 1))
    }
    fclose(fh)
    unlink(fname)
}

void function ECOS_vec_export(real scalar fh, real vector v)
{
    colvector C
    C = bufio()
    fbufput(C, fh, "%4bu", length(v))
    fbufput(C, fh, "%4bu", v)
}

string scalar ECOS_code (real scalar exitcode) {
    if (exitcode == 0 ) return("Optimal solution found")
    if (exitcode == 1 ) return("Certificate of primal infeasibility found")
    if (exitcode == 2 ) return("Certificate of dual infeasibility found")
    if (exitcode == 10) return("Optimal solution found subject to reduced tolerances")
    if (exitcode == 11) return("Certificate of primal infeasibility found subject to reduced tolerances")
    if (exitcode == 12) return("Certificate of dual infeasibility found subject to reduced tolerances")
    if (exitcode == -1) return("Maximum number of iterations reached")
    if (exitcode == -2) return("Numerical problems (unreliable search direction)")
    if (exitcode == -3) return("Numerical problems (slacks or multipliers outside cone)")
    if (exitcode == -4) return("Interrupted by signal or CTRL-C")
    if (exitcode == -7) return("Unknown problem in solver")
}
end
