#ifndef EULER_EXACT_H
#define EULER_EXACT_H

#include"Types.h"

int solve_euler_riemann_exact(
    double rho_L, double u_L, double p_L,
    double rho_R, double u_R, double p_R,
    double gamma, double t, 
    const double *x, int n, double x0,
    double *rho, double *u, double *p);

void io_write_Euler_exact( double *x, double *rho, double *u, double *p,
                          int n, const char *filename);

#endif // EULER_EXACT_H