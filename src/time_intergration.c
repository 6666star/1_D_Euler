#include "time_intergration.h"
#include <string.h>
#include<math.h>

void rk3_init_time(TimeInfo *time,double t_tend,double dt){
    time->current_time = 0.0;
    time->end_time = t_tend;
    time->dt = dt;
}

void rk3_advance_time(DGSolver* solver,FluxFunction flux_func){
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

        if (step_count%10000==0){
            printf("  Step %d, time = %.6f\n", step_count, time->current_time);
        }
        
        spatial_compute_rhs_local(solver, sol->uh, sol->rhs_temp, flux_func);
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < DIM_PK; j++) {
                sol->uh_stage1[i][j] = sol->uh[i][j] + 
                                       (1.0/3.0) * time->dt * sol->rhs_temp[i][j];
            }
        }
        
        spatial_compute_rhs_local(solver, sol->uh_stage1, sol->rhs_temp, flux_func);
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < DIM_PK; j++) {
                sol->uh_stage2[i][j] = sol->uh[i][j] + 
                                       (2.0/3.0) * time->dt * sol->rhs_temp[i][j];
            }
        }
        
        spatial_compute_rhs(solver, sol->uh_stage2, sol->rhs_temp, flux_func);
        spatial_compute_rhs(solver, sol->uh, sol->rhs, flux_func);
        
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < DIM_PK; j++) {
                sol->uh[i][j] = sol->uh[i][j] + 
                               (1.0/4.0) * time->dt * sol->rhs[i][j] + 
                               (3.0/4.0) * time->dt * sol->rhs_temp[i][j];
            }
        }

    }
    
}