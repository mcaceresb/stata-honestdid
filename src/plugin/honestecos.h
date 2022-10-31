#ifndef HONESTECOS
#define HONESTECOS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#define HONESTECOS_VERSION "0.1.1"

#define HONESTECOS_CHAR(cvar, len)                \
    char *(cvar) = malloc(sizeof(char) * (len)); \
    memset ((cvar), '\0', sizeof(char) * (len))

#define HONESTECOS_PWMAX(a, b) ( (a) > (b) ? (a) : (b) )
#define HONESTECOS_PWMIN(a, b) ( (a) > (b) ? (b) : (a) )

idxint     honestecos_read_int(char **bufptr);
ST_retcode honestecos_read_vector(char **bufptr, pfloat **x, idxint *nn);
ST_retcode honestecos_read_intvector(char **bufptr, idxint **x, idxint *nn);
ST_retcode honestecos_read_matrix(char **bufptr, pfloat **A, idxint **i, idxint **p, idxint *A_nnz);
ST_retcode honestecos();

#endif
