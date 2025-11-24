#include "initial_conditions.h"
#include<math.h>
#include<string.h>
#include <stdio.h>


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

void ic_init_solution(DGSolver* solver){
    memset(solver->solution.init, 0, sizeof(solver->solution.init));

    for(int i = 0; i < NX; i++){
        for(int j = 0; j <NUM_GAUSS_POINTS; j++){
            double x = solver->mesh.cell_centers[i] + 
                       solver->mesh.half_cell_size * solver->quad.points[j];
            
            solver->solution.init[i][j][0] = U1(x);
            solver->solution.init[i][j][1] = U2(x);
            solver->solution.init[i][j][2] = U3(x);  // 使用全局gamma
        }
        if(i==NX-1||i==NX-2){
            printf("Initial condition set completed.\n");
               for(int j = 0; j <NUM_GAUSS_POINTS; j++){
                printf("%f  %f  %f\n",solver->solution.init[i][j][0],solver->solution.init[i][j][1],solver->solution.init[i][j][2]);
               }
        }
    }
}

double U1(double x) { return rho_0(x); }
double U2(double x) { return rho_0(x) * u_0(x); }
double U3(double x) { 
    return p_0(x) / (GAMMA - 1.0) + 0.5 * rho_0(x) * u_0(x) * u_0(x);
}
 
void ic_L2_projection(DGSolver* solver) {
    memset(solver->solution.uh, 0, sizeof(solver->solution.uh));

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < DIM_PK; j++) {   
            double int_rho = 0.0;
            double int_rhou = 0.0;
            double int_E = 0.0;
            for (int g = 0; g < NUM_GAUSS_POINTS; g++) {
                
                int_rho+= 0.5*solver->solution.init[i][g][0] * 
                                                solver->basis.phi[g][j] * 
                                                solver->quad.weights[g];
                
                int_rhou+= 0.5*solver->solution.init[i][g][1] * 
                                                solver->basis.phi[g][j] * 
                                                solver->quad.weights[g];
                
               int_E+= 0.5*solver->solution.init[i][g][2] * 
                                                solver->basis.phi[g][j] * 
                                                solver->quad.weights[g];
            }  
                     
            solver->solution.uh[i][j][0] = int_rho / solver->basis.mass_matrix[0][j];
            solver->solution.uh[i][j][1] = int_rhou / solver->basis.mass_matrix[0][j];
            solver->solution.uh[i][j][2] = int_E / solver->basis.mass_matrix[0][j];
        }

    }

}

//---------------定义rho u p函数-----------------
double  rho_0(double x){
    if(x<0.0){
        return 1.0;
    }
    else{
        return 1.0;
    }
}

double  u_0(double x){
    if(x<0.0){
        return 10.0;
    }
    else{
        return 0.0;
    }
}

double  p_0(double x){
    if(x<0.0){
        return 1000.0;
    }
    else{
        return 0.01;
    }
}
//---------------------------------------------
