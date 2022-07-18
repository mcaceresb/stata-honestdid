cap mata mata drop ECOS()
cap mata mata drop ECOS_sol()
cap mata mata drop ECOS_obj()
cap mata mata drop ECOS_setup()
cap mata mata drop ECOS_solve()
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
*     result = ECOS(c, G, h, l, q, e, A, b)
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
*        fname = st_tempfilename()
*        ECOS_setup(fname,
*                   n,
*                   m,
*                   p,
*                   l,
*                   ncones,
*                   q,
*                   e,
*                   G,
*                   A,
*                   c,
*                   h,
*                   b)
*        result = ECOS_solve(fname)
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

real scalar ECOS_obj(real vector c,
                     real matrix G,
                     real vector h,
                     real scalar l,
                     real vector q,
                     real scalar e,
                     real matrix A,
                     real vector b)
{
    struct ECOS_workspace_abridged scalar res
    res = ECOS(c, G, h, l, q, e, A, b)
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
                     real vector b)
{
    struct ECOS_workspace_abridged scalar res
    res = ECOS(c, G, h, l, q, e, A, b)
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
                                           real vector b)
{
    string scalar fname
    real scalar ncones
    fname = st_tempfilename()
    ncones = max(q) > 0? length(q): 0
    ECOS_setup(fname,
               length(c),
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
               b)
    return(ECOS_solve(fname))
}

void ECOS_setup(string scalar fname,
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

struct ECOS_workspace_abridged scalar ECOS_solve(string scalar fname) {
    struct ECOS_workspace_abridged scalar work
    stata(sprintf(`"plugin call honestecos_plugin, `"%s"'"', fname))
    ECOS_cleanup(work, fname)
    return(work)
}

void ECOS_cleanup(struct ECOS_workspace_abridged scalar work, string scalar fname) {
    real scalar fh
    colvector C
    fh = fopen(fname, "r")
    C = bufio()
    work.rc = fbufget(C, fh, "%4bu", 1)
    if ( work.rc == 0 ) {
        work.info_exitcode = fbufget(C, fh, "%4bu", 1)
        work.success       = ((work.info_exitcode == 0) | (work.info_exitcode == 10)) 
        work.info_status   = ECOS_code(work.info_exitcode)
        work.info_obj_val  = fbufget(C, fh, "%8z", 1)
        work.solution_x    = fbufget(C, fh, "%8z", fbufget(C, fh, "%4bu", 1))
        work.solution_y    = fbufget(C, fh, "%8z", fbufget(C, fh, "%4bu", 1))
        work.solution_z    = fbufget(C, fh, "%8z", fbufget(C, fh, "%4bu", 1))
    }
    else {
        work.success       = 0
        work.info_obj_val  = .
        work.info_exitcode = .
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
