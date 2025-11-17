#include "initial_conditions.h"
#include<math.h>
#include<string.h>

void ic_init_mesh(MeshInfo*mesh,double xa,double xb,
                 double bc_left,double bc_right,double t_trend){

mesh->x_start=xa;
mesh->x_end=xb;
mesh->bc_left=bc_left;
mesh->bc_right=bc_right;
mesh->end_time=t_trend;
mesh->cell_size=(xb-xa)/NX;
mesh->half_cell_size=0.5*mesh->cell_size;

    for(int i=0;i<NX;i++){
        mesh->cell_centers[i]=xa+(i+1)*mesh->cell_size-mesh->half_cell_size;
    }
}

void ic_init_(DGSolver* solver){
    memset(solver->solution.uh,0,sizeof(solver->solution.init));

    for(int i=0;i<NX;i++){
        for(int j=0;j<DIM_PK;j++){
     double x = solver->mesh.cell_centers[i] + 
                solver->mesh.half_cell_size * solver->quad.points[j];

                solver->solution.init[i][j][0] = F1(x);
                solver->solution.init[i][j][1] = F2(x);
                solver->solution.init[i][j][2] = F3(x);
        }
    }
}
 
void ic_l2_projection(DGSolver* solver) {
    memset(solver->solution.uh, 0, sizeof(solver->solution.uh));

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < DIM_PK; j++) {   
            for (int g = 0; g < NUM_GAUSS_POINTS; g++) {
                
                solver->solution.uh[i][j][0] += 0.5*solver->solution.init[i][g][0] * 
                                                solver->basis.phi[g][j] * 
                                                solver->quad.weights[g];
                
                solver->solution.uh[i][j][1] += 0.5*solver->solution.init[i][g][1] * 
                                                solver->basis.phi[g][j] * 
                                                solver->quad.weights[g];
                
                solver->solution.uh[i][j][2] += 0.5*solver->solution.init[i][g][2] * 
                                                solver->basis.phi[g][j] * 
                                                solver->quad.weights[g];
            }        
            solver->solution.uh[i][j][0] = solver->solution.uh[i][j][0] / solver->basis.mass_matrix[0][j];
            solver->solution.uh[i][j][1] = solver->solution.uh[i][j][1] / solver->basis.mass_matrix[0][j];
            solver->solution.uh[i][j][2] = solver->solution.uh[i][j][2] / solver->basis.mass_matrix[0][j];
        }
    }
}