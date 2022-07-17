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
    HONESTECOS_CHAR(fname, strlen(argv[0]) + 1);
    strcpy(fname, argv[0]);
    if ( strcmp(fname, "_plugin_check") == 0 ) {
        sf_printf("(note: honestecos_plugin v"HONESTECOS_VERSION" successfully loaded)\n");
    }
    else {
        rc = honestecos(fname);
    }
    return(rc);
}

ST_retcode honestecos(char *fname)
{
    ST_retcode rc = 0;
    ST_double *x = NULL, *y = NULL, *z = NULL, obj;
    uint32_t i;
    FILE *fhandle;

    idxint n, m, p, l, ncones, e, nq, nc, nh, nb, exitflag;
    idxint G_nnz, A_nnz;
    pfloat *G_x = NULL, *A_x = NULL;
    idxint *G_i = NULL, *A_i = NULL;
    idxint *G_p = NULL, *A_p = NULL;
    idxint *q = NULL;
    pfloat *c = NULL, *h = NULL, *b = NULL;
    pwork  *work = NULL;

    // Read data
    // ---------

    fhandle = fopen(fname, "rb");

    if ((rc = honestecos_read_int      (fhandle, &n)))      goto exit;
    if ((rc = honestecos_read_int      (fhandle, &m)))      goto exit;
    if ((rc = honestecos_read_int      (fhandle, &p)))      goto exit;
    if ((rc = honestecos_read_int      (fhandle, &l)))      goto exit;
    if ((rc = honestecos_read_int      (fhandle, &ncones))) goto exit;
    if ((rc = honestecos_read_intvector(fhandle, &q, &nq))) goto exit;
    if ((rc = honestecos_read_int      (fhandle, &e)))      goto exit;

    if ((rc = honestecos_read_matrix(fhandle, &G_x, &G_i, &G_p, &G_nnz))) goto exit;
    if ((rc = honestecos_read_matrix(fhandle, &A_x, &A_i, &A_p, &A_nnz))) goto exit;

    if ((rc = honestecos_read_vector(fhandle, &c, &nc))) goto exit;
    if ((rc = honestecos_read_vector(fhandle, &h, &nh))) goto exit;
    if ((rc = honestecos_read_vector(fhandle, &b, &nb))) goto exit;

    if ((rc = ((x = calloc(HONESTECOS_PWMAX(n, 1), sizeof *x)) == NULL))) goto exit;
    if ((rc = ((y = calloc(HONESTECOS_PWMAX(p, 1), sizeof *y)) == NULL))) goto exit;
    if ((rc = ((z = calloc(HONESTECOS_PWMAX(m, 1), sizeof *z)) == NULL))) goto exit;

    fclose (fhandle);

    // Solve Problem
    // -------------

    work = ECOS_setup(n, m, p, l, ncones, q, e,  
                      G_x, G_p, G_i, // Gpr, Gjc, Gir,  
                      A_x, A_p, A_i, // Apr, Ajc, Air,  
                      c, h, b);

    if (work) {
        work->stgs->verbose = 0; // set to 1 to debug
        exitflag = ECOS_solve(work);
    }
    else {
        rc = 17294;
        goto exit;
    }

    fhandle = fopen(fname, "wb");
    rc = rc | (fwrite((uint32_t *) &rc, sizeof(uint32_t), 1, fhandle) != 1);
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

        rc = rc | (fwrite((uint32_t *) &exitflag, sizeof(uint32_t),  1, fhandle) != 1);
        rc = rc | (fwrite(&obj,                   sizeof(ST_double), 1, fhandle) != 1);
        rc = rc | (fwrite((uint32_t *) &n,        sizeof(uint32_t),  1, fhandle) != 1);
        rc = rc | (fwrite(x,                      sizeof *x,         n, fhandle) != n);
        rc = rc | (fwrite((uint32_t *) &p,        sizeof(uint32_t),  1, fhandle) != 1);
        rc = rc | (fwrite(y,                      sizeof *y,         p, fhandle) != p);
        rc = rc | (fwrite((uint32_t *) &m,        sizeof(uint32_t),  1, fhandle) != 1);
        rc = rc | (fwrite(z,                      sizeof *z,         m, fhandle) != m);
    }
    else {
        rc = 17293;
        goto exit;
    }
    fclose (fhandle);

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

ST_retcode honestecos_read_matrix(FILE *fhandle, pfloat **A_x, idxint **A_i, idxint **A_p, idxint *A_nnz)
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

ST_retcode honestecos_read_vector(FILE *fhandle, pfloat **x, idxint *nn)
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


    *nn = (idxint) n;
    for (j = 0; j < n; j++) {
        (*x)[j] = (pfloat) y[j];
    }

exit:
    if (y) free(y);
    return(rc);
}

ST_retcode honestecos_read_intvector(FILE *fhandle, idxint **x, idxint *nn)
{
    ST_retcode rc = 0;
    uint32_t n, j;
    uint32_t *y = NULL;

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

    *nn = (idxint) n;
    for (j = 0; j < n; j++) {
        (*x)[j] = (idxint) y[j];
    }

exit:
    if (y) free(y);
    return(rc);
}

ST_retcode honestecos_read_int(FILE *fhandle, idxint *nn)
{
    ST_retcode rc = 0;
    uint32_t n, xx;
    xx = fread(&n, sizeof n, 1, fhandle);
    rc = (xx != 1);
    if (rc) goto exit;
    *nn = (idxint) n;
exit:
    return(rc);
}
