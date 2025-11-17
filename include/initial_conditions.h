#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

#include"Types.h"

void ic_init_mesh(MeshInfo*mesh,double xa,double xb,
                 double bc_left,double bc_right,double t_tend);
void ic_set_exact_solution(DGSolver*solver);
void ic_L2_projection(DGSolver*solver);

#endif // INITIAL_CONDITIONS_H