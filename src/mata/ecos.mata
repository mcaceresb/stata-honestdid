cap mata mata drop ECOS()
cap mata mata drop ECOS_get()
cap mata mata drop ECOS_sol()
cap mata mata drop ECOS_obj()
cap mata mata drop ECOS_setup()
cap mata mata drop ECOS_solve()
cap mata mata drop ECOS_read()
cap mata mata drop ECOS_cleanup()
cap mata mata drop ECOS_vec_export()
cap mata mata drop ECOS_code()
cap mata mata drop ECOS_workspace_abridged()

* Second-order cone program
* -------------------------
*
* Mata interface to ECOS. This solves problems of the form
*
*    min  c' x
*    s.t. Ax = b
*         Gx <(K) h
*
* where c is a vector conformable with x, A, G are matrices conformable
* with x, and h is a vector of the same dimension as Gx. In this case,
* <(K) indicates that h - Gx belongs to the cone K, and any number of
* first and second order cones may be specified.
*
* Concretely, say G' = [G1' G2'] and h' = [h1' h2'] with h1-G1x and
* h2-G2x belonging to the first and second-order cones, resp. The
* first-order cone is the positive ortahnt R^n_+, so the constraint is
*
*    G1x <= h1
*
* Now a vector (u_1, u_2) in R x R^{n - 1} belongs to the second-order
* cone if
*
*    ||u_2||_2 <= u_1
*
* with ||.||_2 the L_2 or Eucledean norm. So the constraint is
*
*    L_2(h22 - G22x) <= h21 - G21x
*
* For G2' = [G21' G22'] and h2' = [h21 h22'], noting G21 is a row
* vector, G22 is a matrix, h21 is a scalar, and h22 is a column vector.
* In order to specify multiple first-order cone constraints, we simply
* stack multiple matrices in G. In order to specify multiple second-oder
* cone constraints, we also stack the matrices in G but we need to tell
* ECOS which chunks to treat at sepparate.
*
* With this notation in mind,
*
*     result = ECOS(c, G, h, l, q, e, A, b, 1)
*
* where
*
*     c = vector for objective
*     G = matrix with stacked first and second order cone matrices
*     h = vector with right hand of conic constraint
*     l = Number of first-order cones (linear constraints)
*     q = Vector with the dimensions of the second-order cones
*     e = Number of exponential cones (WARNING: Only e = 0 has been tested)
*     A = Matrix for linear equality constraint
*     b = Vector for equality constraint
*
* result will be a structure of type ECOS_workspace_abridged with components
*
* - rc: Return code for the function (0 means success)
* - info_exitcode: ECOS exit code (0 means solution found, 10 means solution at reduced tol)
* - info_status: Description of ECOS exit code
* - info_obj_val: Objective value
* - solution_x: Minimizer
*
*
* Sparse matrices
* ---------------
*
* If you read the ECOS documentation you will note sparse matrices are
* required.  ECOS converts stata matrices to Compressed Sparse Column
* format internally using the OSQP_csc_matrix structure. If you need
* to create your own sparce matrices using the OSQP_csc_matrix structure,
* you can access the internals as follows:
*
*        bufvar = ECOS_setup(n,
*                            m,
*                            p,
*                            l,
*                            ncones,
*                            q,
*                            e,
*                            G,
*                            A,
*                            c,
*                            h,
*                            b)
*        result = ECOS_solve(bufvar, 1)
*
* Where G, A are already of type OSQP_csc_matrix. See the "Internal
* reference for ECOS_setup" below for details.
*
* Reference
* ---------
*
* See
*
* - https://github.com/embotech/ecos/wiki/Usage-from-C
* - https://github.com/embotech/ecos/wiki/Usage-from-Python
* - Alexander Domahidi's PhD Thesis, "Methods and Tools for Embedded Optimization and Control"
*
* Internal reference for ECOS_setup
* ---------------------------------
*
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
    real vector solution_y
    real vector solution_z
    real scalar success
    real scalar rc
}

real matrix ECOS_get(struct ECOS_workspace_abridged scalar work, string scalar element) {
    if ( element == "info_exitcode" ) return(work.info_exitcode);
    if ( element == "info_obj_val"  ) return(work.info_obj_val);
    if ( element == "solution_x"    ) return(work.solution_x);
    if ( element == "solution_y"    ) return(work.solution_y);
    if ( element == "solution_z"    ) return(work.solution_z);
    if ( element == "success"       ) return(work.success);
    if ( element == "rc"            ) return(work.rc);
}

real scalar ECOS_obj(real vector c,
                     real matrix G,
                     real vector h,
                     real scalar l,
                     real vector q,
                     real scalar e,
                     real matrix A,
                     real vector b,
                     | real scalar cleanup,
                     real scalar verbose)
{
    if ( args() < 9  ) cleanup = 0
    if ( args() < 10 ) verbose = 0
    struct ECOS_workspace_abridged scalar res
    res = ECOS(c, G, h, l, q, e, A, b, cleanup, verbose)
    if ( (res.rc == 0) & ((res.info_exitcode == 0) | (res.info_exitcode == 10)) ) {
        return(res.info_obj_val)
    }
    else {
        errprintf("ECOS failed with status: %g, %s\n", res.info_exitcode, res.info_status)
        _error(res.rc? res.rc: 198)
    }
}

real vector ECOS_sol(real vector c,
                     real matrix G,
                     real vector h,
                     real scalar l,
                     real vector q,
                     real scalar e,
                     real matrix A,
                     real vector b,
                     | real scalar cleanup,
                     real scalar verbose)
{
    if ( args() < 9  ) cleanup = 0
    if ( args() < 10 ) verbose = 0
    struct ECOS_workspace_abridged scalar res
    res = ECOS(c, G, h, l, q, e, A, b, cleanup, verbose)
    if ( (res.rc == 0) & ((res.info_exitcode == 0) | (res.info_exitcode == 10)) ) {
        return(res.solution_x)
    }
    else {
        errprintf("ECOS failed with status: %g, %s\n", res.info_exitcode, res.info_status)
        _error(res.rc? res.rc: 198)
    }
}

struct ECOS_workspace_abridged scalar ECOS(real vector c,
                                           real matrix G,
                                           real vector h,
                                           real scalar l,
                                           real vector q,
                                           real scalar e,
                                           real matrix A,
                                           real vector b,
                                           | real scalar cleanup,
                                           real scalar verbose)
{
    if ( args() < 9  ) cleanup = 0
    if ( args() < 10 ) verbose = 0
    string scalar bufvar
    real scalar ncones
    ncones = max(q) > 0? length(q): 0
    bufvar = ECOS_setup(length(c),
                        rows(G),
                        rows(A),
                        l,
                        ncones,
                        q,
                        e,
                        OSQP_csc_convert(G),
                        OSQP_csc_convert(A),
                        c,
                        h,
                        b,
                        verbose)
    return(ECOS_solve(bufvar, cleanup))
}

string scalar ECOS_setup(real scalar n,
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
                         real vector b,
                         | real scalar verbose)
{
    if ( args() < 13 ) verbose = 0

    string scalar buf, bufvar
    real scalar off

    st_numscalar("__honestecos_rc",    .)
    st_numscalar("__honestecos_exit",  .)
    st_numscalar("__honestecos_obj",   .)
    st_numscalar("__honestecos_maxit", length(st_numscalar("__honestecos_maxit"))? st_numscalar("__honestecos_maxit"): .)

    st_matrix("__honestecos_x", J(1, n, .))
    st_matrix("__honestecos_y", J(1, p, .))
    st_matrix("__honestecos_z", J(1, m, .))

    off = 0
    buf = char(0)

    OSQP_int_export(buf, verbose, off)
    OSQP_int_export(buf, n,       off)
    OSQP_int_export(buf, m,       off)
    OSQP_int_export(buf, p,       off)
    OSQP_int_export(buf, l,       off)
    OSQP_int_export(buf, ncones,  off)
    ECOS_vec_export(buf, q,       off)
    OSQP_int_export(buf, e,       off)
    OSQP_csc_export(buf, G,       off)
    OSQP_csc_export(buf, A,       off)
    OSQP_vec_export(buf, c,       off)
    OSQP_vec_export(buf, h,       off)
    OSQP_vec_export(buf, b,       off)

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

struct ECOS_workspace_abridged scalar ECOS_solve(string scalar bufvar, | real scalar cleanup) {
    if ( args() < 2 ) cleanup = 0
    struct ECOS_workspace_abridged scalar work
    stata(sprintf(`"plugin call honestecos_plugin %s, _plugin_run"', bufvar))
    ECOS_read(work, cleanup)
    return(work)
}

void ECOS_read(struct ECOS_workspace_abridged scalar work, | real scalar cleanup) {
    work.rc = st_numscalar("__honestecos_rc")
    if ( work.rc == 0 ) {
        work.info_exitcode = st_numscalar("__honestecos_exit")
        work.success       = ((work.info_exitcode == 0) | (work.info_exitcode == 10)) 
        work.info_status   = ECOS_code(work.info_exitcode)
        work.info_obj_val  = st_numscalar("__honestecos_obj")
        work.solution_x    = st_matrix("__honestecos_x")
        work.solution_y    = st_matrix("__honestecos_y")
        work.solution_z    = st_matrix("__honestecos_z")
    }
    else {
        work.success       = 0
        work.info_obj_val  = .
        work.info_exitcode = .
    }
    if ( cleanup ) ECOS_cleanup()
}

void ECOS_cleanup() {
    real scalar ix
    ix = _st_varindex(st_global("__HONESTBUFVAR"))
    if ( !missing(ix) ) (void) st_dropvar(ix)
    st_global("__HONESTBUFVAR", "")

    st_numscalar("__honestecos_rc",    J(0, 0, .))
    st_numscalar("__honestecos_exit",  J(0, 0, .))
    st_numscalar("__honestecos_obj",   J(0, 0, .))
    st_numscalar("__honestecos_maxit", J(0, 0, .))

    st_matrix("__honestecos_x", J(0, 0, .))
    st_matrix("__honestecos_y", J(0, 0, .))
    st_matrix("__honestecos_z", J(0, 0, .))
}

void function ECOS_vec_export(string scalar buf, real vector v, real scalar off)
{
    colvector C
    C = bufio()
    buf = buf + (4 + 4 * length(v)) * char(0);
    bufput(C, buf, off, "%4bu", length(v)); off = off + 4;
    bufput(C, buf, off, "%4bu", v);         off = off + 4 * length(v);
}

string scalar ECOS_code (real scalar exitcode) {
    string scalar decode
    decode = "Unknown exit code"
    if (exitcode == 0 ) decode = "Optimal solution found"
    if (exitcode == 1 ) decode = "Certificate of primal infeasibility found"
    if (exitcode == 2 ) decode = "Certificate of dual infeasibility found"
    if (exitcode == 10) decode = "Optimal solution found subject to reduced tolerances"
    if (exitcode == 11) decode = "Certificate of primal infeasibility found subject to reduced tolerances"
    if (exitcode == 12) decode = "Certificate of dual infeasibility found subject to reduced tolerances"
    if (exitcode == -1) decode = "Maximum number of iterations reached"
    if (exitcode == -2) decode = "Numerical problems (unreliable search direction)"
    if (exitcode == -3) decode = "Numerical problems (slacks or multipliers outside cone)"
    if (exitcode == -4) decode = "Interrupted by signal or CTRL-C"
    if (exitcode == -7) decode = "Unknown problem in solver"
    return(decode)
}
end
