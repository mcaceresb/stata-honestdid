// gcc -I/home/mauricio/bulk/lib/osqp/include honestdid.c /home/mauricio/bulk/lib/osqp/build/out/libosqp.a -lm -o honestdid && ./honestdid

#include "osqp.h"

int main(int argc, char **argv) {
    int i;

    // Load problem data
    c_float P_x[10] = {1.0632808e-03, -5.7471800e-05,  7.3610400e-04, -2.5761580e-04, -3.1207080e-04,  8.1683460e-04, -1.6259000e-05,  6.3565800e-05, -1.9438180e-04,  6.1791900e-04,};
    c_int P_nnz = 10;
    c_int P_i[10] = {4, 4, 5, 4, 5, 6, 4, 5, 6, 7,};
    c_int P_p[9] = {0, 0, 0, 0, 0, 1, 3, 6, 10, };
    c_float q[8] = {0, 0, 0, 0, -1.051178e-04, 7.224761e-05, 1.285712e-04, 2.395732e-04,};

    // Matrix A
    c_int A_nnz     = 32;
    c_float A_x[32] = {-1., -1., -1., -1., -1., -1., -1., -1.,  1.,  1.,  1.,  1.,  1., -1., -1., -1., -1.,  1.,  1.,  1.,  1., -1., -1., -1.,  1.,  1.,  1., -1., -1.,  1.,  1., -1.,};
    c_int A_i[32]   = {1, 5, 2, 6, 3, 7, 4, 8, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 2, 3, 4, 6, 7, 8, 0, 3, 4, 7, 8, 0, 4, 8,};
    c_int A_p[9]    = { 0,  2,  4,  6,  8, 17, 24, 29, 32};
    c_float u[9]    = {1, 0, 0, 0, 0, 0, 0, 0, 0,};
    c_float l[9]    = {1, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY,};
    c_int n         = 8;
    c_int m         = 9;

    // Exitflag
    c_int exitflag = 0;

    // Workspace structures
    OSQPWorkspace *work;
    OSQPSettings  *settings = (OSQPSettings *) c_malloc(sizeof(OSQPSettings));
    OSQPData      *data     = (OSQPData *)     c_malloc(sizeof(OSQPData));

    // Populate data
    if (data) {
        data->n = n;
        data->m = m;
        data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
        data->q = q;
        data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
        data->l = l;
        data->u = u;
    }

    // Define solver settings as default
    if (settings) {
        osqp_set_default_settings(settings);
    }

    // Setup workspace
    exitflag = osqp_setup(&work, data, settings);
    printf("setup: %d\n", exitflag);

    // Solve Problem
    exitflag = osqp_solve(work);
    printf("solution: %d\n", exitflag);

printf("little waffles: %.6f\n", work->info->obj_val);
printf("debug: so why the segfault????\n");
printf("debug: \t\t%.9g\n", work->solution->x[0]);
printf("debug: \t\t%.9g\n", work->solution->x[1]);
printf("debug: \t\t%.9g\n", work->solution->x[2]);
printf("debug: \t\t%.9g\n", work->solution->x[3]);
printf("debug: \t\t%.9g\n", work->solution->x[4]);
printf("debug: \t\t%.9g\n", work->solution->x[5]);
printf("debug: \t\t%.9g\n", work->solution->x[6]);
printf("debug: \t\t%.9g\n", work->solution->x[7]);

    // Cleanup
    osqp_cleanup(work);
    if (data) {
        if (data->A) c_free(data->A);
        if (data->P) c_free(data->P);
        c_free(data);
    }
    if (settings) c_free(settings);


    return exitflag;
}
