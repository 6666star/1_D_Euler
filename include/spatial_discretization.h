// ==================== spatial_discretization.h ====================
#ifndef SPATIAL_DISCRETIZATION_H
#define SPATIAL_DISCRETIZATION_H

#include "Types.h"

// 标准DG离散
void spatial_compute_rhs(DGSolver* solver,
                        const  double uh[NX][DIM_PK][3],
                        double rhs[NX][DIM_PK][3]);
// 局部DG离散
void spatial_compute_rhs_local(DGSolver* solver,
                        const  double uh[NX][DIM_PK][3],
                        double rhs[NX][DIM_PK][3],
                        FluxFunction flux_func);

// 通量函数
double HLL_flux(double uL,double uR,double fL,double fR,double SR,double SL);

double LF_flux(double uL,double uR,double fL,double fR,double SR,double SL);

void wavespeed(double u1, double u2, double u3, DGSolver* solver, double s[2]);

double func_2(double u1,double u2,double u3,DGSolver* solver);

double func_1(double u1,double u2,double u3,DGSolver* solver);

double func_0(double u1,double u2,double u3,DGSolver* solver);


#endif //SPATIAL_DISCRETIZATION_H