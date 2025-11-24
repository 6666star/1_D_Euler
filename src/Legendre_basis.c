#include"Legendre_basis.h"
#include<string.h>

void basis_init(BasisFunctions* basis, const GaussQuadrature* quad,double hx){
  memset(basis,0,sizeof(BasisFunctions));

  double hx_inV=2.0/hx;       //网格1/2*hx的倒数

  for(int i=0;i<NUM_GAUSS_POINTS;i++){
    double xi=quad->points[i];

    basis->phi[i][0]=1.0;
    basis->phi[i][1]=xi;
    basis->phi[i][2]=xi*xi-(1.0/3.0);
    //basis->phi[i][3]=xi*xi*xi-(3.0/5.0)*xi;

    basis->phi_x[i][0]=0.0;
    basis->phi_x[i][1]=hx_inV;
    basis->phi_x[i][2]=2.0*xi*hx_inV;
    //basis->phi_x[i][3]=3.0*xi*xi*hx_inV-(3.0/5.0)*hx_inV;
  }

    basis->phi_right[0][0] = 1.0;
    basis->phi_right[0][1] = 1.0;
    basis->phi_right[0][2] = 2.0 / 3.0;
    //basis->phi_right[0][3] = 2.0 / 5.0;

    basis->phi_left[0][0] = 1.0;
    basis->phi_left[0][1] = -1.0;
    basis->phi_left[0][2] = 2.0 / 3.0;
    //basis->phi_left[0][3] = -2.0 / 5.0;

    basis->mass_matrix[0][0] = 1.0;
    basis->mass_matrix[0][1] = 1.0 / 3.0;
    basis->mass_matrix[0][2] = 4.0 / 45.0;
    //basis->mass_matrix[0][3] = 4.0 / 175.0;

    
    //    printf("Mass matrix diagonal:\n");
    // for(int i=0;i<NUM_GAUSS_POINTS;i++){
    //   printf("%f  ",basis->phi_x[i][0]);
    // }
}

double basis_eval_legendre(int order,double xi){
    switch (order)
    {
    case 0: return 1.0;
    case 1: return xi;
    case 2: return xi*xi-(1.0/3.0);
    //case 3: return xi*xi*xi-(3.0/5.0)*xi;
    default:return 0.0;
    }
}

double basisi_eval_legendre_derivative(int order,double xi){
    switch (order) {
    case 0: return 0.0;
    case 1: return 1.0;
    case 2: return 2.0*xi;
    //case 3: return 3.0*xi*xi-(3.0/5.0);
    default:return 0.0;
    }
}