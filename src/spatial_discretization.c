#include "spatial_discretization.h"
#include <string.h>
#include <math.h>
#include <stdio.h>

//---------------Euler方程通量函数-----------------
double func_0(double u1,double u2,double u3,DGSolver* solver){
    return u2;
}

double func_1(double u1,double u2,double u3,DGSolver* solver){
   double gamma = GAMMA;

  double p = (gamma - 1.0) * (u3 - 0.5 * (u2 * u2) / u1);
  return (u2 * u2) / u1 + p;
}

double func_2(double u1,double u2,double u3,DGSolver* solver){
  double gamma = GAMMA;

  double p = (gamma - 1.0) * (u3 - 0.5 * (u2 * u2) / u1);
  return u2 / u1 * (u3 + p);
}
//-------------------------------------------------

//---------------计算特征速度---------------------
void wavespeed(double u1, double u2, double u3, DGSolver* solver, double s[2])
{
    double gamma = GAMMA;

    double velocity = u2 / u1;
    double pressure = (gamma - 1.0) * (u3 - 0.5 * (u2 * u2) / u1);
    double sound_speed = sqrt(gamma * pressure / u1);

    s[0] = velocity - sound_speed;  // 左波速
    s[1] = velocity + sound_speed;  // 右波速
}
//-------------------------------------------------


//---------------通量函数---------------------
//---------------Lax-Friedrichs---------------------
double LF_flux(double uL,double uR,double fL,double fR,double SR,double SL){
    double alpha=fmax(fabs(SR),fabs(SL));

    return 0.5*(fL+fR)-0.5*alpha*(uR-uL);
}
//-------------------------------------------------
//---------------HLL通量---------------------------
double HLL_flux(double uL,double uR,double fL,double fR,double SR,double SL){
    if(SL>=0){
        return fL;
    }else if(SR<=0){
        return fR;
    }else{
        return (SR*fL-SL*fR+SR*SL*(uR-uL))/(SR-SL);
    }
}
//-------------------------------------------------


void spatial_compute_rhs(DGSolver* solver,
                        const  double uh[NX][DIM_PK][3],
                        double rhs[NX][DIM_PK][3]){

    FluxArrays* flux=&solver->flux_work;                
    BasisFunctions* basis=&solver->basis;
    MeshInfo* mesh=&solver->mesh;
    GaussQuadrature* quad=&solver->quad;

    memset(rhs,0,sizeof(double)*NX*DIM_PK*3);
    memset(flux->uh_gauss,0,sizeof(flux->uh_gauss));
    memset(flux->uh_left,0,sizeof(flux->uh_left));
    memset(flux->uh_right,0,sizeof(flux->uh_right));
    memset(flux->flux,0,sizeof(flux->flux));

    //------------------   周期边界 -------------------
    for(int i=0;i<NX;i++){
        for(int j=0;j<DIM_PK;j++){
            flux->uh_boundary[i+1][j][0]=uh[i][j][0];
            flux->uh_boundary[i+1][j][1]=uh[i][j][1];
            flux->uh_boundary[i+1][j][2]=uh[i][j][2];
        }
    }

       if (mesh->bc_left == 1.0 && mesh->bc_right == 1.0) {
        for (int j = 0; j < DIM_PK; j++) {
            flux->uh_boundary[0][j][0] = uh[NX-1][j][0];
            flux->uh_boundary[NX+1][j][0] = uh[0][j][0];
            flux->uh_boundary[0][j][1] = uh[NX-1][j][1];
            flux->uh_boundary[NX+1][j][1] = uh[0][j][1];
            flux->uh_boundary[0][j][2] = uh[NX-1][j][2];
            flux->uh_boundary[NX+1][j][2] = uh[0][j][2];
            }
        }

        if (mesh->bc_left == 2.0 && mesh->bc_right == 2.0) {
        for (int j = 0; j < DIM_PK; j++) {
            flux->uh_boundary[0][j][0] = uh[0][j][0];
            flux->uh_boundary[NX+1][j][0] = uh[NX-1][j][0];
            flux->uh_boundary[0][j][1] = uh[0][j][1];
            flux->uh_boundary[NX+1][j][1] = uh[NX-1][j][1];
            flux->uh_boundary[0][j][2] = uh[0][j][2];
            flux->uh_boundary[NX+1][j][2] = uh[NX-1][j][2];
            }
        }
    //-------------------------------------------------

    //------------------- uh积分投影到guass点上 -------------------
    for(int i=0;i<NX;i++){
        for(int j=0;j<DIM_PK;j++){
            for(int k=0;k<NUM_GAUSS_POINTS;k++){
                flux->uh_gauss[i][k][0]+=uh[i][j][0]*basis->phi[k][j];
                flux->uh_gauss[i][k][1]+=uh[i][j][1]*basis->phi[k][j];
                flux->uh_gauss[i][k][2]+=uh[i][j][2]*basis->phi[k][j];
            }
        }
    }

    //------------------------------------------------------
//    for(int i=0;i<3;i++){
//     printf("位置%d:\n",i);
//     for(int j=0;j<NX;j++){
//         for(int k=0;k<NUM_GAUSS_POINTS;k++){
//             printf("%f ",flux->uh_gauss[j][k][i]);
//         }
//      }
//    }
 
    //-----------------计算空间积分------------------
    for(int i=0;i<NX;i++){
        for(int j=0;j<DIM_PK;j++){
            for(int k=0;k<NUM_GAUSS_POINTS;k++){
                rhs[i][j][0]+=0.5*quad->weights[k]*func_0(flux->uh_gauss[i][k][0],
                                           flux->uh_gauss[i][k][1],
                                           flux->uh_gauss[i][k][2],solver)
                                           *basis->phi_x[k][j];

                rhs[i][j][1]+=0.5*quad->weights[k]*func_1(flux->uh_gauss[i][k][0],
                                           flux->uh_gauss[i][k][1],
                                           flux->uh_gauss[i][k][2],solver)
                                           *basis->phi_x[k][j];

                rhs[i][j][2]+=0.5*quad->weights[k]*func_2(flux->uh_gauss[i][k][0],
                                           flux->uh_gauss[i][k][1],
                                           flux->uh_gauss[i][k][2],solver)
                                           *basis->phi_x[k][j];
            }
        }
    }
    //------------------------------------------------

//    for(int i=0;i<3;i++){
//     printf("模态位置%d:\n",i);
//     for(int j=0;j<NX;j++){
//         for(int k=0;k<DIM_PK;k++){
//             printf("%15f ",rhs[j][k][i]);
//         }
//         printf("\n");
//     }
//    }
    //-------------计算数值通量--------------------
    for(int i=0;i<NX+1;i++){
        for(int j=0;j<DIM_PK;j++){
           flux->uh_left[i][0][0]+=flux->uh_boundary[i+1][j][0]*basis->phi_left[0][j];
           flux->uh_left[i][0][1]+=flux->uh_boundary[i+1][j][1]*basis->phi_left[0][j];
           flux->uh_left[i][0][2]+=flux->uh_boundary[i+1][j][2]*basis->phi_left[0][j];

           flux->uh_right[i][0][0]+=flux->uh_boundary[i][j][0]*basis->phi_right[0][j];
           flux->uh_right[i][0][1]+=flux->uh_boundary[i][j][1]*basis->phi_right[0][j];
           flux->uh_right[i][0][2]+=flux->uh_boundary[i][j][2]*basis->phi_right[0][j];
        }
    }

//        for(int i=0;i<3;i++){
//     printf("uh_right and uh_left%d:\n",i);
//     for(int j=0;j<NX+1;j++){
        
//             printf("%15f  %15f\n",flux->uh_left[j][0][i],flux->uh_right[j][0][i]);
    
//     }
//    }

    for(int i=0;i<NX+1;i++){
        double f_left[3],f_right[3];
 
        f_right[0]=func_0(flux->uh_left[i][0][0],flux->uh_left[i][0][1],flux->uh_left[i][0][2],solver);
        f_right[1]=func_1(flux->uh_left[i][0][0],flux->uh_left[i][0][1],flux->uh_left[i][0][2],solver);
        f_right[2]=func_2(flux->uh_left[i][0][0],flux->uh_left[i][0][1],flux->uh_left[i][0][2],solver);

        f_left[0]=func_0(flux->uh_right[i][0][0],flux->uh_right[i][0][1],flux->uh_right[i][0][2],solver); 
        f_left[1]=func_1(flux->uh_right[i][0][0],flux->uh_right[i][0][1],flux->uh_right[i][0][2],solver);
        f_left[2]=func_2(flux->uh_right[i][0][0],flux->uh_right[i][0][1],flux->uh_right[i][0][2],solver);
        

       //   printf("界面%d左右状态变量%d: %f %f %f\n",i, f_left[0],f_left[1],f_left[2]);
        
       double SL[2];
       double SR[2];

       wavespeed(flux->uh_left[i][0][0],flux->uh_left[i][0][1],flux->uh_left[i][0][2],solver,SL);
       wavespeed(flux->uh_right[i][0][0],flux->uh_right[i][0][1],flux->uh_right[i][0][2],solver,SR);

       double S_min = (SL[0] < SR[0]) ? SL[0] : SR[0];
       double S_max = (SL[1] > SR[1]) ? SL[1] : SR[1];
        
      // printf("界面%d左波速%f右波速%f\n",i,S_min,S_max);
        for(int j=0;j<3;j++){
            flux->flux[i][0][j]=LF_flux(flux->uh_right[i][0][j],flux->uh_left[i][0][j],
                                     f_right[j],f_left[j],S_max,S_min);                                
        }
           
    }


//    for(int i=0;i<3;i++){
//     printf("fhat计算%d:\n",i);
//     for(int j=0;j<NX+1;j++){
//             printf("%15f\n ",flux->flux[j][0][i]);
    
//     }
//    }


    //-----------------------------------------------   

    //----------------------组装各项--------------------
    for(int i=0;i<NX;i++){
        for(int j=0;j<DIM_PK;j++){
            for(int k=0;k<3;k++){
                rhs[i][j][k]-=(1.0/mesh->cell_size)*(basis->phi_right[0][j]*flux->flux[i+1][0][k]
                              -basis->phi_left[0][j]*flux->flux[i][0][k]);
            }
        }

        for(int j=0;j<DIM_PK;j++){
            rhs[i][j][0]=rhs[i][j][0]/basis->mass_matrix[0][j];
            rhs[i][j][1]=rhs[i][j][1]/basis->mass_matrix[0][j];
            rhs[i][j][2]=rhs[i][j][2]/basis->mass_matrix[0][j];
        }
    }


//        for(int i=0;i<3;i++){
//     printf("模态位置%d:\n",i);
//     for(int j=0;j<NX;j++){
//         for(int k=0;k<DIM_PK;k++){
//             printf("%15f ",rhs[j][k][i]);
//         }
//         printf("\n");
//     }
//    }
}
    



// void spatial_compute_rhs_local(DGSolver* solver,
//                         const  double uh[NX][DIM_PK][3],
//                         double rhs[NX][DIM_PK][3],
//                         FluxFunction flux_func){

//     FluxArrays* flux=&solver->flux_work;                
//     BasisFunctions* basis=&solver->basis;
//     MeshInfo* mesh=&solver->mesh;
//     GaussQuadrature* quad=&solver->quad;

//     memset(rhs,0,sizeof(double)*NX*DIM_PK*3);
//     memset(flux->uh_gauss,0,sizeof(flux->uh_gauss));
//     memset(flux->uh_left,sizeof(flux->uh_left));
//     memset(flux->uh_right,sizeof(flux->uh_right));
//     memset(flux->flux,0,sizeof(flux->flux));

//     //------------------   周期边界 -------------------
//     for(int i=0;i<NX;i++){
//         for(int j=0;j<NUM_GAUSS_POINTS;j++){
//             flux->uh_boundary[i+1][j][0]=uh[i][j][0];
//             flux->uh_boundary[i+1][j][0]=uh[i][j][1];
//             flux->uh_boundary[i+1][j][0]=uh[i][j][2];
//         }
//     }

//        if (mesh->bc_left == 1.0 && mesh->bc_right == 1.0) {
//         for (int j = 0; j < DIM_PK; j++) {
//             flux->uh_boundary[0][j][0] = uh[NX-1][j][0];
//             flux->uh_boundary[NX+1][j][0] = uh[0][j][0];
//             flux->uh_boundary[0][j][1] = uh[NX-1][j][1];
//             flux->uh_boundary[NX+1][j][1] = uh[0][j][1];
//             flux->uh_boundary[0][j][2] = uh[NX-1][j][2];
//             flux->uh_boundary[NX+1][j][2] = uh[0][j][2];
//         }
//     }
//     //-------------------------------------------------

//     //------------------- uh积分投影到guass点上 -------------------
//     for(int i=0;i<NX;i++){
//         for(int j=0;j<DIM_PK;j++){
//             for(int k=0;k<NUM_GAUSS_POINTS;k++){
//                 flux->uh_gauss[i][k][0]+=uh[i][j][0]*basis->phi[k][j];
//                 flux->uh_gauss[i][k][1]+=uh[i][j][1]*basis->phi[k][j];
//                 flux->uh_gauss[i][k][2]+=uh[i][j][2]*basis->phi[k][j];
//             }
//         }
//     }
//     //------------------------------------------------------
 
//     //-----------------计算空间积分------------------
//     for(int i=0;i<NX;i++){
//         for(int j=0;j<DIM_PK;j++){
//             for(int k=0;k<NUM_GAUSS_POINTS;k++){
//                 rhs[i][j][0]+=0.5*quad->weights[k]*func_0(uh[i][j][0],uh[i][j][1],uh[i][j][2],solver)
//                               *basis->phi_x[k][j];
//                 rhs[i][j][1]+=0.5*quad->weights[k]*func_1(uh[i][j][0],uh[i][j][1],uh[i][j][2],solver)
//                               *basis->phi_x[k][j];
//                 rhs[i][j][2]+=0.5*quad->weights[k]*func_2(uh[i][j][0],uh[i][j][1],uh[i][j][2],solver)
//                               *basis->phi_x[k][j];
//             }
//         }
//     }
//     //------------------------------------------------

//     //-------------计算数值通量--------------------
//     for(int i=0;i<=NX+1;i++){
//         for(int j=0;j<DIM_PK;j++){
//            flux->uh_left[i][0][0]+=flux->uh_boundary[i][0][0]*basis->phi_left[0][j];
//            flux->uh_left[i][0][1]+=flux->uh_boundary[i][0][1]*basis->phi_left[0][j];
//            flux->uh_left[i][0][2]+=flux->uh_boundary[i][0][2]*basis->phi_left[0][j];

//            flux->uh_right[i][0][0]+=flux->uh_boundary[i][0][0]*basis->phi_right[0][j];
//            flux->uh_right[i][0][1]+=flux->uh_boundary[i][0][1]*basis->phi_right[0][j];
//            flux->uh_right[i][0][2]+=flux->uh_boundary[i][0][2]*basis->phi_right[0][j];
//         }
//     }

//     for(i=0;i<Nx+1;i++){
//         double f_left[3],f_right[3];
//         double
 
//         f_left[0]=func_0(flux->uh_left[i][0][0],flux->uh_left[i][0][1],flux->uh_left[i][0][2],solver);
//         f_left[1]=func_1(flux->uh_left[i][0][0],flux->uh_left[i][0][1],flux->uh_left[i][0][2],solver);
//         f_left[2]=func_2(flux->uh_left[i][0][0],flux->uh_left[i][0][1],flux->uh_left[i][0][2],solver);

//         f_right[0]=func_0(flux->uh_right[i][0][0],flux->uh_right[i][0][1],flux->uh_right[i][0][2],solver); 
//         f_right[1]=func_1(flux->uh_right[i][0][0],flux->uh_right[i][0][1],flux->uh_right[i][0][2],solver);
//         f_right[2]=func_2(flux->uh_right[i][0][0],flux->uh_right[i][0][1],flux->uh_right[i][0][2],solver);
    

//         flux->flux[i][0][0]=func_1(flux->uh_left[i][0][0],flux->uh_left[i][0][1],flux->uh_left[i][0][2],solver);
//         flux->flux[i][0][1]=func_2(flux->uh_left[i][0][0],flux->uh_left[i][0][1],flux->uh_left[i][0][2],solver);
//         flux->flux[i][0][2]=func_3(flux->uh_left[i][0][0],flux->uh_left[i][0][1],flux->uh_left[i][0][2],solver);

//         flux->flux[i][1][0]=func_1(flux->uh_right[i][0][0],flux->uh_right[i][0][1],flux->uh_right[i][0][2],solver);
//         flux->flux[i][1][1]=func_2(flux->uh_right[i][0][0],flux->uh_right[i][0][1],flux->uh_right[i][0][2],solver);
//         flux->flux[i][1][2]=func_3(flux->uh_right[i][0][0],flux->uh_right[i][0][1],flux->uh_right[i][0][2],solver);
        
//     }

//     //-----------------------------------------------   

//     //----------------------组装各项--------------------
//     for(int i=0;i<NX;i++){
//         for(int j=0;j<DIM_PK;j++){
//             for(int k=0;k<3;k++){
//                 rhs[i][j][k]-=(1.0/mesh->cell_size)*basis->phi_right[0][j]*flux->flux[i+1][0][k]
//                               -basis->phi_left[0][j]*flux->flux[i][1][k];
//             }
//         }

//         for(int j=0;j<DIM_PK;j++){
//             rhs[i][j][0]=rhs[i][j][0]/basis->mass_matrix[0][j];
//             rhs[i][j][1]=rhs[i][j][1]/basis->mass_matrix[0][j];
//             rhs[i][j][2]=rhs[i][j][2]/basis->mass_matrix[0][j];
//         }
//     }
//    //--------------------------------------------- 
// }  

