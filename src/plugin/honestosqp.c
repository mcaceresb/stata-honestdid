// 17290 -      OOM init
// 17291 - OSQP OOM data
// 17292 - OSQP OOM settings
// 17293 - OSQP OOM work
// 17294 - ECOS OOM work

#include "stplugin.h"
#include "osqp.h"
#include "honestosqp.h"
#include "sf_printf.c"

STDLL stata_call(int argc, char * argv[])
{
    ST_retcode rc = 0;
    HONESTOSQP_CHAR(todo, strlen(argv[0]) + 1);
    strcpy(todo, argv[0]);
    if ( strcmp(todo, "_plugin_check") == 0 ) {
        sf_printf("(note: honestosqp_plugin v"HONESTOSQP_VERSION" successfully loaded)\n");
    }
    else if ( strcmp(todo, "_plugin_run") == 0 ) {
        rc = honestosqp();
    }
    else {
        sf_printf("unknown option %s\n", todo);
        rc = 198;
    }
    return(rc);
}

ST_retcode honestosqp()
{
    ST_retcode rc = 0;
    char *buf = NULL, *bufptr;
    ST_double *x = NULL, obj;
    ST_int bytes = SF_sdatalen(1, 1);

    c_int verbose, n, m, i;
    c_int P_nnz, A_nnz;
    c_float *P_x = NULL, *A_x = NULL;
    c_int   *P_i = NULL, *A_i = NULL;
    c_int   *P_p = NULL, *A_p = NULL;
    c_float *q = NULL, *u = NULL, *l = NULL;

    // Workspace structures
    // --------------------

    OSQPWorkspace *work = NULL;
    OSQPSettings  *settings = (OSQPSettings *) c_malloc(sizeof(OSQPSettings));
    OSQPData      *data     = (OSQPData *)     c_malloc(sizeof(OSQPData));
    data->A = NULL;
    data->P = NULL;

    // Load problem data
    // -----------------

    buf = malloc(bytes);
    if ( buf == NULL ) {
        rc = 17290;
        goto exit;
    }
    bufptr = buf;

    if ((rc = SF_strldata(1, 1, buf, bytes)) != bytes) {
        goto exit;
    }
    else rc = 0;
    verbose = honestosqp_read_int(&bufptr);

    if ((rc = honestosqp_read_matrix(&bufptr, &P_x, &P_i, &P_p, &P_nnz))) goto exit;
    if ((rc = honestosqp_read_matrix(&bufptr, &A_x, &A_i, &A_p, &A_nnz))) goto exit;

    if ((rc = honestosqp_read_vector(&bufptr, &q, &n))) goto exit;
    if ((rc = honestosqp_read_vector(&bufptr, &u, &m))) goto exit;
    if ((rc = honestosqp_read_vector(&bufptr, &l, &m))) goto exit;

    if ((rc = ((x = calloc(n, sizeof *x)) == NULL))) goto exit;

    for (i = 0; i < m; i++) {
        if ( SF_is_missing(l[i]) ) l[i] = -OSQP_INFTY;
        if ( SF_is_missing(u[i]) ) u[i] = OSQP_INFTY;
    }

    // Debugging
    // ---------

// printf("bytes = %d\n", bytes);
// printf("n = %lld\n", n);
// for(i = 0; i < n; i++) {
//     printf("\tq[%lld] = %.9f\n", i, q[i]);
// }
//     printf("\n\n");
// printf("m = %lld\n", m);
// for(i = 0; i < m; i++) {
//     if ( l[i] == -OSQP_INFTY && u[i] == OSQP_INFTY ) {
//         printf("\t(l, u)[%lld] = -infty, infty\n", i);
//     }
//     else if ( l[i] == -OSQP_INFTY ) {
//         printf("\t(l, u)[%lld] = -infty, %.9f\n", i, u[i]);
//     }
//     else if ( u[i] == OSQP_INFTY ) {
//         printf("\t(l, u)[%lld] = %.9f, infty\n", i, l[i]);
//     }
//     else {
//         printf("\t(l, u)[%lld] = %.9f, %.9f\n", i, l[i], u[i]);
//     }
// }
//     printf("\n\n");
// printf("P_nnz = %lld\n", P_nnz);
// for(i = 0; i < P_nnz; i++) {
//     printf("\tP[%lld] = %.9f\n", i, P_x[i]);
// }
//     printf("\n\n");
// printf("A_nnz = %lld\n", A_nnz);
// for(i = 0; i < A_nnz; i++) {
//     printf("\tP[%lld] = %.9f\n", i, A_x[i]);
// }
//     printf("\n\n");

    // Populate data
    // -------------

    if (data) {
        data->n = n;
        data->m = m;
        data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
        data->q = q;
        data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
        data->l = l;
        data->u = u;
    }
    else {
        rc = 17291;
        goto exit;
    }

    // Define solver settings as default
    // ---------------------------------

    if (settings) {
        osqp_set_default_settings(settings);
        // Default parameters to match CVXR:
        // https://cvxr.rbind.io/cvxr_examples/cvxr_solver-parameters/
        settings->verbose      = verbose;
        settings->eps_abs      = 1e-5;
        settings->eps_rel      = 1e-5;
        settings->eps_prim_inf = 1e-4;
        settings->eps_dual_inf = 1e-4;
        settings->max_iter     = 10000;

        // rho is adapted every 25 iterations (by default) if some
        // fraction of elapsed time is greater than the setup time.  I
        // suppose the thought is that it's not necessary to adapt rho
        // when the algorithm is converging quickly; however, this leads
        // to non-deterministic output in certain cases when rho is not
        // adapted (this is the culprit of GH issue #7, and the reason
        // why it's hard to replicate).
        settings->adaptive_rho_interval = 25;
    }
    else {
        rc = 17292;
        goto exit;
    }

    // Solve Problem
    // -------------

    if ( data && settings ) {
        if ((rc = osqp_setup(&work, data, settings))) goto exit;
        if ((rc = osqp_solve(work))) goto exit;
        for (i = 0; i < n; i++) {
            x[i] = work->solution->x[i];
        }
        obj = work->info->obj_val;
    }
    else {
        rc = 17291;
        goto exit;
    }

    rc = rc | (SF_scal_save("__honestosqp_rc", (ST_double) rc));
    if (work) {
        rc = rc | (SF_macro_save("___honestosqp_status", work->info->status));
        rc = rc | (SF_scal_save("__honestosqp_obj", obj));
        for (i = 0; i < n; i++) {
            rc = rc | SF_mat_store("__honestosqp_x", 1, i + 1, x[i]);
        }
    }
    else {
        rc = 17293;
        goto exit;
    }

    // Cleanup
    // -------

exit:

    if (work) osqp_cleanup(work);
    if (data) {
        if (data->A) c_free(data->A);
        if (data->P) c_free(data->P);
        c_free(data);
    }
    if (settings) c_free(settings);

    if (P_x) free(P_x);
    if (P_i) free(P_i);
    if (P_p) free(P_p);

    if (A_x) free(A_x);
    if (A_i) free(A_i);
    if (A_p) free(A_p);

    if (x)   free(x);
    if (q)   free(q);
    if (u)   free(u);
    if (l)   free(l);

    if (buf) free(buf);

    return(rc);
}

ST_retcode honestosqp_read_matrix(char **bufptr, c_float **A_x, c_int **A_i, c_int **A_p, c_int *A_nnz)
{
    ST_retcode rc = 0;
    uint32_t nnz, ni, np, j;
    uint32_t *i = NULL, *p = NULL;
    ST_double *A = NULL;

    memmove(&nnz, *bufptr, sizeof nnz); *bufptr += (sizeof nnz);
    memmove(&ni,  *bufptr, sizeof ni);  *bufptr += (sizeof ni);
    memmove(&np,  *bufptr, sizeof np);  *bufptr += (sizeof np);

    rc = rc | (((*A_x) = calloc(nnz, sizeof **A_x)) == NULL);
    rc = rc | (((*A_i) = calloc(ni,  sizeof **A_i)) == NULL);
    rc = rc | (((*A_p) = calloc(np,  sizeof **A_p)) == NULL);
    rc = rc | ((A      = calloc(nnz, sizeof *A))    == NULL);
    rc = rc | ((i      = calloc(ni,  sizeof *i))    == NULL);
    rc = rc | ((p      = calloc(np,  sizeof *p))    == NULL);
    if (rc) {
        rc = 17290;
        goto exit;
    }

    memmove(A, *bufptr, (sizeof *A) * nnz); *bufptr += (sizeof *A) *  nnz;
    memmove(i, *bufptr, (sizeof *i) * ni);  *bufptr += (sizeof *i) *  ni;
    memmove(p, *bufptr, (sizeof *p) * np);  *bufptr += (sizeof *p) *  np;

    *A_nnz = (c_int) nnz;
    for (j = 0; j < nnz; j++) {
        (*A_x)[j] = (c_float) A[j];
    }

    for (j = 0; j < ni; j++) {
        (*A_i)[j] = (c_int) i[j];
    }

    for (j = 0; j < np; j++) {
        (*A_p)[j] = (c_int) p[j];
    }

exit:
    if (A) free(A);
    if (i) free(i);
    if (p) free(p);
    return(rc);
}

ST_retcode honestosqp_read_vector(char **bufptr, c_float **x, c_int *nn)
{
    ST_retcode rc = 0;
    uint32_t n, j;
    ST_double *y = NULL;

    memmove(&n, *bufptr, sizeof n); *bufptr += (sizeof n);

    rc = (((*x) = calloc(n, sizeof **x)) == NULL);
    rc = ((  y  = calloc(n, sizeof  *y)) == NULL);
    if (rc) {
        rc = 17290;
        goto exit;
    }

    memmove(y, *bufptr, (sizeof *y) * n); *bufptr += (sizeof *y) * n;
    *nn = (c_int) n;
    for (j = 0; j < n; j++) {
        (*x)[j] = (c_float) y[j];
    }

exit:
    if (y) free(y);
    return(rc);
}

c_int honestosqp_read_int(char **bufptr)
{
    uint32_t n;
    memmove(&n, *bufptr, sizeof n); *bufptr += (sizeof n); 
    return((c_int) n);
}
