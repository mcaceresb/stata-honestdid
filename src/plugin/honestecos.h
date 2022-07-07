#ifndef HONESTECOS
#define HONESTECOS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#define HONESTECOS_VERSION "0.1.0"

#define HONESTECOS_CHAR(cvar, len)                \
    char *(cvar) = malloc(sizeof(char) * (len)); \
    memset ((cvar), '\0', sizeof(char) * (len))

#define HONESTECOS_PWMAX(a, b) ( (a) > (b) ? (a) : (b) )
#define HONESTECOS_PWMIN(a, b) ( (a) > (b) ? (b) : (a) )

ST_retcode honestecos_read_int(FILE *fhandle, idxint *nn);
ST_retcode honestecos_read_vector(FILE *fhandle, pfloat **x, idxint *nn);
ST_retcode honestecos_read_intvector(FILE *fhandle, idxint **x, idxint *nn);
ST_retcode honestecos_read_matrix(FILE *fhandle, pfloat **A, idxint **i, idxint **p, idxint *A_nnz);
ST_retcode honestecos(char *fname);

#endif
