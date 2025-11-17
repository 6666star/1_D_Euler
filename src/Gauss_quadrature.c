#include"Gauss_quadrature.h"
#include<string.h>

void guass_init_quadrature(GaussQuadrature* quad){
    quad->points[0] = -0.9061798459386639927976269;
    quad->points[1] = -0.5384693101056830910363144;
    quad->points[2] = 0.0;
    quad->points[3] = 0.5384693101056830910363144;
    quad->points[4] = 0.9061798459386639927976269;

    quad->weights[0] = 0.2369268850561890875142640;
    quad->weights[1] = 0.4786286704993664680412915;
    quad->weights[2] = 0.5688888888888888888888889;
    quad->weights[3] = 0.4786286704993664680412915;
    quad->weights[4] = 0.2369268850561890875142640;
}