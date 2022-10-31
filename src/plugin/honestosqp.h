#ifndef HONESTOSQP
#define HONESTOSQP

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#define HONESTOSQP_VERSION "0.1.0"

#define HONESTOSQP_CHAR(cvar, len)                \
    char *(cvar) = malloc(sizeof(char) * (len)); \
    memset ((cvar), '\0', sizeof(char) * (len))

#define HONESTOSQP_PWMAX(a, b) ( (a) > (b) ? (a) : (b) )
#define HONESTOSQP_PWMIN(a, b) ( (a) > (b) ? (b) : (a) )

ST_retcode honestosqp_read_vector(char **bufptr, c_float **x, c_int *nn);
ST_retcode honestosqp_read_matrix(char **bufptr, c_float **A, c_int **i, c_int **p, c_int *A_nnz);
c_int      honestosqp_read_int(char **bufptr);
ST_retcode honestosqp();

#endif
