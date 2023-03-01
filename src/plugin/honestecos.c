// 17290 -      OOM init
// 17291 - OSQP OOM data
// 17292 - OSQP OOM settings
// 17293 - OSQP OOM work
// 17294 - ECOS OOM work

#include "stplugin.h"
#include "ecos.h"
#include "honestecos.h"
#include "sf_printf.c"

STDLL stata_call(int argc, char * argv[])
{
    ST_retcode rc = 0;
    HONESTECOS_CHAR(todo, strlen(argv[0]) + 1);
    strcpy(todo, argv[0]);
    if ( strcmp(todo, "_plugin_check") == 0 ) {
        sf_printf("(note: honestecos_plugin v"HONESTECOS_VERSION" successfully loaded)\n");
    }
    else if ( strcmp(todo, "_plugin_run") == 0 ) {
        rc = honestecos();
    }
    else {
        sf_printf("unknown option %s\n", todo);
        rc = 198;
    }
    return(rc);
}

ST_retcode honestecos()
{
    ST_retcode rc = 0;
    char *buf = NULL, *bufptr;
    ST_double *x = NULL, *y = NULL, *z = NULL, obj, maxit;
    uint32_t i;
    ST_int bytes = SF_sdatalen(1, 1);

    idxint verbose, n, m, p, l, ncones, e, nq, nc, nh, nb, exitflag;
    idxint G_nnz, A_nnz;
    pfloat *G_x = NULL, *A_x = NULL;
    idxint *G_i = NULL, *A_i = NULL;
    idxint *G_p = NULL, *A_p = NULL;
    idxint *q = NULL;
    pfloat *c = NULL, *h = NULL, *b = NULL;
    pwork  *work = NULL;

    // Read data
    // ---------

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

    verbose = honestecos_read_int(&bufptr);
    n       = honestecos_read_int(&bufptr);
    m       = honestecos_read_int(&bufptr);
    p       = honestecos_read_int(&bufptr);
    l       = honestecos_read_int(&bufptr);
    ncones  = honestecos_read_int(&bufptr);
    if ((rc = honestecos_read_intvector(&bufptr, &q, &nq))) goto exit;
    e       = honestecos_read_int(&bufptr);
 

    if ((rc = SF_scal_use("__honestecos_maxit", &maxit))) goto exit;

    if ((rc = honestecos_read_matrix(&bufptr, &G_x, &G_i, &G_p, &G_nnz))) goto exit;
    if ((rc = honestecos_read_matrix(&bufptr, &A_x, &A_i, &A_p, &A_nnz))) goto exit;

    if ((rc = honestecos_read_vector(&bufptr, &c, &nc))) goto exit;
    if ((rc = honestecos_read_vector(&bufptr, &h, &nh))) goto exit;
    if ((rc = honestecos_read_vector(&bufptr, &b, &nb))) goto exit;

    if ((rc = ((x = calloc(HONESTECOS_PWMAX(n, 1), sizeof *x)) == NULL))) goto exit;
    if ((rc = ((y = calloc(HONESTECOS_PWMAX(p, 1), sizeof *y)) == NULL))) goto exit;
    if ((rc = ((z = calloc(HONESTECOS_PWMAX(m, 1), sizeof *z)) == NULL))) goto exit;

    // Solve Problem
    // -------------

    work = ECOS_setup(n, m, p, l, ncones, q, e,  
                      G_x, G_p, G_i, // Gpr, Gjc, Gir,  
                      A_x, A_p, A_i, // Apr, Ajc, Air,  
                      c, h, b);

    if (work) {
        // Default parameters to match CVXR:
        // https://cvxr.rbind.io/cvxr_examples/cvxr_solver-parameters/
        work->stgs->verbose = verbose;
        work->stgs->reltol  = 1e-8;
        work->stgs->abstol  = 1e-8;
        work->stgs->feastol = 1e-8;
        work->stgs->maxit   = (idxint) (SF_is_missing(maxit)? 100: maxit);
        exitflag = ECOS_solve(work);
    }
    else {
        rc = 17294;
        goto exit;
    }

    rc = rc | (SF_scal_save("__honestecos_rc", (ST_double) rc));
    rc = rc | (SF_scal_save("__honestecos_exit", (ST_double) exitflag));
    if (work) {
        // Proofing against corner instance where n, p, m = 0 somehow
        x[0] = SV_missval;
        y[0] = SV_missval;
        z[0] = SV_missval;
        obj  = SV_missval;

        if ( (exitflag != 0) && (exitflag != 10)  ) {
            for (i = 0; i < n; i++) {
                x[i] = SV_missval;
            }
            for (i = 0; i < p; i++) {
                y[i] = SV_missval;
            }
            for (i = 0; i < m; i++) {
                z[i] = SV_missval;
            }
        }
        else {
            obj = 0;
            for (i = 0; i < n; i++) {
                x[i] = work->x[i];
                obj += c[i] * x[i];
            }
            for (i = 0; i < p; i++) {
                y[i] = work->y[i];
            }
            for (i = 0; i < m; i++) {
                z[i] = work->z[i];
            }
        }

        // Proofing against corner instance where n, p, m = 0 somehow
        n = HONESTECOS_PWMAX(n, 1);
        p = HONESTECOS_PWMAX(p, 1);
        m = HONESTECOS_PWMAX(m, 1);

        rc = rc | (SF_scal_save("__honestecos_obj", obj));
        for (i = 0; i < n; i++) {
            rc = rc | SF_mat_store("__honestecos_x", 1, i + 1, x[i]);
        }
        for (i = 0; i < p; i++) {
            rc = rc | SF_mat_store("__honestecos_y", 1, i + 1, y[i]);
        }
        for (i = 0; i < m; i++) {
            rc = rc | SF_mat_store("__honestecos_z", 1, i + 1, z[i]);
        }
    }
    else {
        rc = 17293;
        goto exit;
    }

    // Cleanup
    // -------

exit:

    if (work) {
        ECOS_cleanup(work, 0);
    }

    if (G_x) free(G_x);
    if (G_i) free(G_i);
    if (G_p) free(G_p);

    if (A_x) free(A_x);
    if (A_i) free(A_i);
    if (A_p) free(A_p);

    if (q) free(q);
    if (c) free(c);
    if (h) free(h);
    if (b) free(b);

    if (x) free(x);
    if (y) free(y);
    if (z) free(z);

    return(rc);
}

ST_retcode honestecos_read_matrix(char **bufptr, pfloat **A_x, idxint **A_i, idxint **A_p, idxint *A_nnz)
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

    *A_nnz = (idxint) nnz;
    for (j = 0; j < nnz; j++) {
        (*A_x)[j] = (pfloat) A[j];
    }

    for (j = 0; j < ni; j++) {
        (*A_i)[j] = (idxint) i[j];
    }

    for (j = 0; j < np; j++) {
        (*A_p)[j] = (idxint) p[j];
    }

exit:
    if (A) free(A);
    if (i) free(i);
    if (p) free(p);
    return(rc);
}

ST_retcode honestecos_read_vector(char **bufptr, pfloat **x, idxint *nn)
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
    *nn = (idxint) n;
    for (j = 0; j < n; j++) {
        (*x)[j] = (pfloat) y[j];
    }

exit:
    if (y) free(y);
    return(rc);
}

ST_retcode honestecos_read_intvector(char **bufptr, idxint **x, idxint *nn)
{
    ST_retcode rc = 0;
    uint32_t n, j;
    uint32_t *y = NULL;

    memmove(&n, *bufptr, sizeof n); *bufptr += (sizeof n);

    rc = (((*x) = calloc(n, sizeof **x)) == NULL);
    rc = ((  y  = calloc(n, sizeof  *y)) == NULL);
    if (rc) {
        rc = 17290;
        goto exit;
    }

    memmove(y, *bufptr, (sizeof *y) * n); *bufptr += (sizeof *y) * n;
    *nn = (idxint) n;
    for (j = 0; j < n; j++) {
        (*x)[j] = (idxint) y[j];
    }

exit:
    if (y) free(y);
    return(rc);
}

idxint honestecos_read_int(char **bufptr)
{
    uint32_t n;
    memmove(&n, *bufptr, sizeof n); *bufptr += (sizeof n); 
    return((idxint) n);
}
