#ifndef LEGENDERE_BASIS_H
#define LEGENDERE_BASIS_H

#include"Types.h"
void basis_init(BasisFunctions*basis,const GaussQuadrature* quad ,double hx);
double basis_eval_legendre(int order,double xi);
double basisi_eval_legendre_derivative(int order,double xi);

#endif // LEGENDERE_BASIS_H