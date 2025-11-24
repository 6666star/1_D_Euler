#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

#include"Types.h"

void ic_init_mesh(MeshInfo*mesh,double xa,double xb,
                 double bc_left,double bc_right,double t_tend);

void ic_L2_projection(DGSolver*solver);

void ic_init_solution(DGSolver* solver);

double rho_0(double x);
double u_0(double x);
double p_0(double x);

double U1(double x);
double U2(double x);
double U3(double x);

#endif // INITIAL_CONDITIONS_H