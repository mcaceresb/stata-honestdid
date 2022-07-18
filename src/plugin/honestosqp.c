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
    HONESTOSQP_CHAR(fname, strlen(argv[0]) + 1);
    strcpy(fname, argv[0]);
    if ( strcmp(fname, "_plugin_check") == 0 ) {
        sf_printf("(note: honestosqp_plugin v"HONESTOSQP_VERSION" successfully loaded)\n");
    }
    else {
        rc = honestosqp(fname);
    }
    return(rc);
}

ST_retcode honestosqp(char *fname)
{
    ST_retcode rc = 0;
    FILE *fhandle;
    ST_double *x = NULL, obj;

    c_int n, m, i;
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

    // Load problem data
    // -----------------

    fhandle = fopen(fname, "rb");

    if ((rc = honestosqp_read_matrix(fhandle, &P_x, &P_i, &P_p, &P_nnz))) goto exit;
    if ((rc = honestosqp_read_matrix(fhandle, &A_x, &A_i, &A_p, &A_nnz))) goto exit;

    if ((rc = honestosqp_read_vector(fhandle, &q, &n))) goto exit;
    if ((rc = honestosqp_read_vector(fhandle, &u, &m))) goto exit;
    if ((rc = honestosqp_read_vector(fhandle, &l, &m))) goto exit;

    if ((rc = ((x = calloc(n, sizeof *x)) == NULL))) goto exit;

    for (i = 0; i < m; i++) {
        if ( SF_is_missing(l[i]) ) l[i] = -OSQP_INFTY;
        if ( SF_is_missing(u[i]) ) u[i] = OSQP_INFTY;
    }

    fclose (fhandle);

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
        settings->verbose      = 0;
        settings->eps_abs      = 1e-5;
        settings->eps_rel      = 1e-5;
        settings->eps_prim_inf = 1e-4;
        settings->eps_dual_inf = 1e-4;
        settings->max_iter     = 10000;
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

    fhandle = fopen(fname, "wb");

    rc = rc | (fwrite((uint32_t *) &rc, sizeof(uint32_t), 1, fhandle) != 1);
    if (work) {
        rc = rc | (fwrite(work->info->status, 32,                1, fhandle) != 1);
        rc = rc | (fwrite(&obj,               sizeof(ST_double), 1, fhandle) != 1);
        rc = rc | (fwrite((uint32_t *) &n,    sizeof(uint32_t),  1, fhandle) != 1);
        rc = rc | (fwrite(x,                  sizeof *x,         n, fhandle) != n);
    }
    else {
        rc = 17293;
        goto exit;
    }

    fclose (fhandle);

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

    if (fname) free(fname);

    return(rc);
}

ST_retcode honestosqp_read_matrix(FILE *fhandle, c_float **A_x, c_int **A_i, c_int **A_p, c_int *A_nnz)
{
    ST_retcode rc = 0;
    uint32_t nnz, ni, np, j;
    uint32_t *i = NULL, *p = NULL;
    ST_double *A = NULL;

    rc = rc | (fread(&nnz, sizeof nnz, 1, fhandle) != 1);
    rc = rc | (fread(&ni,  sizeof ni,  1, fhandle) != 1);
    rc = rc | (fread(&np,  sizeof np,  1, fhandle) != 1);
    if (rc) goto exit;

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

    rc = rc | (fread(A, sizeof *A, nnz, fhandle) != nnz);
    rc = rc | (fread(i, sizeof *i, ni,  fhandle) != ni);
    rc = rc | (fread(p, sizeof *p, np,  fhandle) != np);
    if (rc) goto exit;

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

ST_retcode honestosqp_read_vector(FILE *fhandle, c_float **x, c_int *nn)
{
    ST_retcode rc = 0;
    uint32_t n, j;
    ST_double *y = NULL;

    rc = (fread(&n, sizeof n, 1, fhandle) != 1);
    if (rc) goto exit;

    rc = (((*x) = calloc(n, sizeof **x)) == NULL);
    rc = ((  y  = calloc(n, sizeof  *y)) == NULL);
    if (rc) {
        rc = 17290;
        goto exit;
    }

    rc = (fread(y, sizeof *y, n, fhandle) != n);
    if (rc) goto exit;


    *nn = (c_int) n;
    for (j = 0; j < n; j++) {
        (*x)[j] = (c_float) y[j];
    }

exit:
    if (y) free(y);
    return(rc);
}
