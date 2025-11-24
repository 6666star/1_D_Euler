#include "time_integration.h"
#include "spatial_discretization.h"
#include "TVD_Limiter.h"
#include <string.h>
#include<math.h>
#include <stdio.h>

void rk3_init_time(TimeInfo *time,double t_tend,double dt){
    time->current_time = 0.0;
    time->end_time = t_tend;
    time->dt = dt;
}

void rk3_advance_time(DGSolver* solver){
    TimeInfo *time = &solver->time;
    SolutionArrays *sol = &solver->solution;
    
    int step_count=0;
    memset(sol->rhs,0,sizeof(sol->rhs));
    
    printf("Starting time integration start .....\n");

    while (time->current_time<time->end_time){
        if(time->current_time+time->dt>time->end_time){ 
            time->dt=time->end_time-time->current_time;
            time->current_time=time->end_time;
        }else{
            time->current_time+=time->dt;
            step_count++;
        }

        // if(step_count==10){
        //     break;
        // }

        if (step_count%10000==0){
            printf("  Step %d, time = %.6f\n", step_count, time->current_time);
        }
        
        spatial_compute_rhs(solver, sol->uh, sol->rhs_temp);
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < DIM_PK; j++) {
                for(int k=0;k<3;k++){
                sol->uh_stage1[i][j][k] = sol->uh[i][j][k] + 
                                         time->dt * sol->rhs_temp[i][j][k];
                }
            }
        }


       TVD_Limiter(solver, sol->uh_stage1);


        
        spatial_compute_rhs(solver, sol->uh_stage1, sol->rhs_temp);
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < DIM_PK; j++) {
                for(int k=0;k<3;k++){
                sol->uh_stage2[i][j][k] =  (3.0/4.0)*sol->uh[i][j][k] +
                                           (1.0/4.0)* sol->uh_stage1[i][j][k]+
                                           (1.0/4.0) * time->dt* sol->rhs_temp[i][j][k];
                }
            }
        }

       TVD_Limiter(solver, sol->uh_stage2);
        
       spatial_compute_rhs(solver, sol->uh_stage2, sol->rhs_temp);
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < DIM_PK; j++) {
                for(int k=0;k<3;k++){
                sol->uh[i][j][k] = (1.0/3.0) *sol->uh[i][j][k] + 
                                   (2.0/3.0) * sol->uh_stage2[i][j][k] + 
                                   (2.0/3.0) * time->dt * sol->rhs_temp[i][j][k];
                }
            }
        }

         TVD_Limiter(solver, sol->uh);

        double min_val = sol->uh[0][0][0]; 
        double max_val = sol->uh[0][0][0]; 
        for (int i = 1; i < NX; ++i) {
            if (sol->uh[i][0][0] < min_val) {
                min_val = sol->uh[i][1][1];
            }
            if (sol->uh[i][0][0] > max_val) {
                max_val = sol->uh[i][1][1];
            }
      }
       

    printf("时间=%fmin_val=%fmax_evl=%f %f \n",time->current_time,min_val,max_val);

//     for(int i=0;i<3;i++){
//     printf("模态位置%d:\n",i);
//     for(int j=0;j<NX;j++){
//         for(int k=0;k<DIM_PK;k++){
//             printf("%15f ",sol->uh[j][k][i]);
//         }
//         printf("\n");
//     }
//    }


    }
    
}