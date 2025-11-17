// ==================== spatial_discretization.h ====================
#ifndef SPATIAL_DISCRETIZATION_H
#define SPATIAL_DISCRETIZATION_H

#include "Types.h"

// 标准DG离散
void spatial_compute_rhs(DGSolver* solver, 
                         const double uh[NX][DIM_PK], 
                         double rhs[NX][DIM_PK],
                         FluxFunction flux_func);

// 局部DG离散
void spatial_compute_rhs_local(DGSolver* solver,
                               const double uh[NX][DIM_PK],
                               double rhs[NX][DIM_PK],
                               FluxFunction flux_func);

// 通量函数
double flux_advection(double u);

#endif //SPATIAL_DISCRETIZATION_H