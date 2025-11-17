#include <stdio.h>
#include <stdlib.h>
#include "dg_types.h"
#include "gauss_quadrature.h"
#include "basis_functions.h"
#include "initial_conditions.h"
#include "time_integration.h"
#include "spatial_discretization.h"
#include "io_utils.h"

int main(int argc, char** argv) {
    (void)argc;  // 避免未使用参数警告
    (void)argv;
    
    // 开始计时
    double start_time = io_get_wall_time();
    
    // 初始化求解器
    DGSolver solver = {0};
    
    printf("=== DG Solver for 1D Euler systems Equation ===\n");
    printf("Grid points: %d\n", NX);
    printf("Polynomial order: %d\n", POLY_ORDER);
    printf("Gauss quadrature points: %d\n\n", NUM_GAUSS_POINTS);
    
    // 1. 初始化高斯积分
    printf("Initializing Gauss quadrature...\n");
    gauss_init_quadrature(&solver.quad);
    
    // 2. 初始化网格
    printf("Setting up mesh...\n");
    ic_init_mesh(&solver.mesh, 0.0, 2.0 * PI, 1.0, 1.0, 2.0 * PI);
    
    // 3. 初始化基函数
    printf("Computing basis functions...\n");
    basis_init(&solver.basis, &solver.quad, solver.mesh.cell_size);
    
    // 4. 设置初始条件
    printf("Setting initial conditions...\n");
    ic_set_exact_solution(&solver);
    ic_l2_projection(&solver);
    
    // 5. 时间步进
    printf("Starting time integration...\n");
    double dt = CFL_NUMBER * solver.mesh.cell_size;
    rk3_init_time(&solver.time, solver.mesh.end_time, dt);
    
    rk3_advance(&solver, flux_advection);
    
    // 6. 输出结果
    printf("\nWriting solution to file...\n");
    io_write_solution(&solver, "DG_solution.dat");
    
    // 7. 计算误差
    printf("Computing L2 error...\n");
    io_compute_error(&solver);
    
    // 结束计时
    double elapsed = io_get_wall_time() - start_time;
    printf("\nWall time elapsed: %.6f seconds\n", elapsed);
    
    printf("\n=== Simulation completed successfully ===\n");
    
    return 0;
}