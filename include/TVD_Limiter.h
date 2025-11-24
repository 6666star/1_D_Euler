//===============TVD_Limter.h=============
#ifndef TVD_LIMTIER_H
#define TVD_LIMTIER_H

#include"Types.h"

void TVD_Limiter(DGSolver* solver,double uh[NX][DIM_PK][3]);

void minmod(double a[3], double b[3], double c[3], double result[3], 
                   double hx, double m);

void mat_inv_3x3(double mat[3][3], double inv[3][3]);

void mat_vec_mult(double mat[3][3], double vec[3], double result[3]) ;

#endif // TVD_LIMTIER_H
