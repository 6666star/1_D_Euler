//==================== apps/main.c (修改版) ====================
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Types.h"
#include "Gauss_quadrature.h"
#include "Legendre_basis.h"
#include "initial_conditions.h"
#include "spatial_discretization.h"
#include "time_integration.h"
#include "TVD_Limiter.h"
#include "io_utils.h"
#include "Euler_exact.h"

// 辅助函数：检查数值有效性
int check_valid_solution(DGSolver* solver) {
    int has_error = 0;
    
    for (int i = 0; i < NX; i++) {
        double rho = solver->solution.uh[i][0][0];
        double rho_u = solver->solution.uh[i][0][1];
        double E = solver->solution.uh[i][0][2];
        
        if (isnan(rho) || isnan(rho_u) || isnan(E)) {
            printf("  ✗ NaN detected at cell %d\n", i);
            has_error = 1;
        }
        if (isinf(rho) || isinf(rho_u) || isinf(E)) {
            printf("  ✗ Inf detected at cell %d\n", i);
            has_error = 1;
        }
        if (rho < 0) {
            printf("  ✗ Negative density at cell %d: rho = %.6e\n", i, rho);
            has_error = 1;
        }
        
        double u = rho_u / rho;
        double p = (GAMMA - 1.0) * (E - 0.5 * rho * u * u);
        if (p < 0) {
            printf("  ✗ Negative pressure at cell %d: p = %.6e\n", i, p);
            has_error = 1;
        }
    }
    
    return !has_error;
}

int main(int argc, char** argv) {
    (void)argc;
    (void)argv;
    
    // 开始计时
    double start_time = io_get_wall_time();
    
    printf("\n");
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║    DG Solver for 1D Euler Equations - Init Test          ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n");
    printf("\n");
    
    // 初始化求解器
    DGSolver solver = {0};
    
    // ==========================================================
    // 打印配置信息
    // ==========================================================
    printf("Configuration:\n");
    printf("  Grid points (NX):            %d\n", NX);
    printf("  Polynomial order:            %d\n", POLY_ORDER);
    printf("  DOF per cell (DIM_PK):       %d\n", DIM_PK);
    printf("  Gauss quadrature points:     %d\n", NUM_GAUSS_POINTS);
    printf("  CFL number:                  %.3f\n", CFL_NUMBER);
    printf("  Adiabatic index (gamma):     %.6f\n", GAMMA);
    printf("\n");
    
    // ==========================================================
    // 1. 初始化高斯积分
    // ==========================================================
    printf("[1/5] Initializing Gauss quadrature...\n");
    gauss_init_quadrature(&solver.quad);
    
    // 验证高斯积分
    double sum_weights = 0.0;
    for (int i = 0; i < NUM_GAUSS_POINTS; i++) {
        sum_weights += solver.quad.weights[i];
    }
    printf("      Sum of weights = %.15f (expected: 2.0)\n", sum_weights);
    if (fabs(sum_weights - 2.0) < 1e-12) {
        printf("      ✓ Gauss quadrature initialized correctly\n");
    } else {
        printf("      ✗ ERROR: Gauss quadrature weights incorrect!\n");
        return 1;
    }
    printf("\n");
    
    // ==========================================================
    // 2. 初始化网格
    // ==========================================================
    printf("[2/5] Setting up mesh...\n");
    double xa = -1.0;
    double xb = 1.0;
    double t_end = 0.01;  // Sod激波管终止时间
  
    ic_init_mesh(&solver.mesh, xa, xb, 2.0, 2.0, t_end);
    
    printf("      Domain: [%.3f, %.3f]\n", solver.mesh.x_start, solver.mesh.x_end);
    printf("      Cell size (dx): %.10f\n", solver.mesh.cell_size);
    printf("      End time: %.3f\n", solver.mesh.end_time);
    printf("      ✓ Mesh initialized\n");
    printf("\n");
    
    // ==========================================================
    // 3. 初始化基函数
    // ==========================================================
    printf("[3/5] Computing basis functions...\n");
    basis_init(&solver.basis, &solver.quad, solver.mesh.cell_size);

    printf("      Number of basis functions: %d\n", DIM_PK);
    printf("      Mass matrix diagonal:\n");
    // for (int j = 0; j < DIM_PK && j < 3; j++) {
    //     printf("        M[%d] = %.10f\n", j, solver.basis.mass_matrix[j]);
    // }
    if (DIM_PK > 3) printf("        ...\n");
    printf("      ✓ Basis functions initialized\n");
    printf("\n");
    
    // ==========================================================
    // 4. 设置初始条件
    // ==========================================================
    printf("[4/5] Setting initial conditions (Sod shock tube)...\n");
    printf("      Left state  (x < 0.5): rho=1.0, u=0.0, p=1.0\n");
    printf("      Right state (x > 0.5): rho=0.125, u=0.0, p=0.1\n");
    
    // 调用你写的初始化函数
    ic_init_solution(&solver);
    
    printf("      ✓ Initial values set at quadrature points\n");
    printf("\n");
    
    // ==========================================================
    // 5. L2投影
    // ==========================================================
    printf("[5/5] Performing L2 projection...\n");
    
    // 调用你写的L2投影函数
    ic_L2_projection(&solver);
    
    printf("      ✓ L2 projection completed\n");
    printf("\n");

    // ==========================================================
    // 5.时间推进
    // ==========================================================
      printf("[6/6] RK3 time intgeral...\n");
    double alpha=1.0;
    double S[2];
    for(int i=0;i<NX;i++){
            wavespeed(solver.solution.uh[i][0][0]
                     ,solver.solution.uh[i][0][1]
                     ,solver.solution.uh[i][0][2],&solver,S); 
             double S_max = (S[0] > S[1]) ? S[0] : S[1];
           if (S_max>alpha)
           {
            alpha=S_max;
           }
    }

    double dt = CFL_NUMBER * solver.mesh.cell_size/alpha;
    rk3_init_time(&solver.time, solver.mesh.end_time, dt);

    rk3_advance_time(&solver);
    // ==========================================================
    // 验证结果
    // ==========================================================
    printf("Verification:\n");
    printf("  Checking solution validity...\n");
    if (check_valid_solution(&solver)) {
        printf("  ✓ All values are physical (no NaN/Inf, rho>0, p>0)\n");
    } else {
        printf("  ✗ ERROR: Invalid values detected!\n");
        return 1;
    }
    
    // 检查守恒量
    double total_mass = 0.0;
    double total_momentum = 0.0;
    double total_energy = 0.0;
    
    for (int i = 0; i < NX; i++) {
        double dx = solver.mesh.cell_size;
        total_mass += solver.solution.uh[i][0][0] * dx;
        total_momentum += solver.solution.uh[i][0][1] * dx;
        total_energy += solver.solution.uh[i][0][2] * dx;
    }
    
    printf("\n  Conservation check:\n");
    printf("    Total mass:     %.10f\n", total_mass);
    printf("    Total momentum: %.10e (should be ~0)\n", total_momentum);
    printf("    Total energy:   %.10f\n", total_energy);
    
    // 理论值
    double L = xb - xa;
    double rho_L = 1.0, p_L = 1.0;
    double rho_R = 0.125, p_R = 0.1;
    double E_L = p_L / (GAMMA - 1.0);
    double E_R = p_R / (GAMMA - 1.0);
    double theory_mass = 0.5 * L * (rho_L + rho_R);
    double theory_energy = 0.5 * L * (E_L + E_R);
    
    printf("    Expected mass:   %.10f (error: %.2e)\n", 
           theory_mass, fabs(total_mass - theory_mass));
    printf("    Expected energy: %.10f (error: %.2e)\n", 
           theory_energy, fabs(total_energy - theory_energy));
    
    // ==========================================================
    // 输出结果
    // ==========================================================
    printf("\nWriting output files...\n");
    
    // 创建输出目录
    system("mkdir -p output");
    
    // 使用你的io_utils.c中的函数
    io_write_solution(&solver, "output/euler_initial.dat");
    io_write_modal_coefficients(&solver, "output/euler_modal.dat");
    io_write_high_resolution(&solver, "output/euler_highres.dat");
    
    printf("  ✓ output/euler_initial.dat\n");
    printf("  ✓ output/euler_modal.dat\n");
    printf("  ✓ output/euler_highres.dat\n");
    printf("\n");
    // ==========================================================
    // 计算解析解并输出
    // ==========================================================
    double rho_exactL = 1.0, u_exactL = 10.0, p_exactL = 1000.0;
    double rho_exactR = 1.0, u_exactR = 0.0,  p_exactR = 0.01;
    double t = 0.01;
    double x0 = 0.0;
    
    // 网格设置
    int n = 1000;
    
    // 分配内存
    double *x = malloc(n * sizeof(double));
    double *rho_exact = malloc(n * sizeof(double));
    double *u_exact = malloc(n * sizeof(double));
    double *p_exact = malloc(n * sizeof(double));
    
    // 生成网格
    for (int i = 0; i < n; i++) {
        x[i] = xa + (xb - xa) * i / (n - 1);
    }

    int status = solve_euler_riemann_exact(
                 rho_exactL, u_exactL, p_exactL,
                 rho_exactR,u_exactR, p_exactR,
                 GAMMA, t, x, n, x0, 
                 rho_exact, u_exact, p_exact);

    if (status == 0) {
       io_write_Euler_exact(x,rho_exact, u_exact, p_exact, n, "output/euler_Exact.dat");
    } else {
       printf("Solver failed!\n");
    }
        // 释放内存
    free(x);
    free(rho_exact);
    free(u_exact);
    free(p_exact);
    // ==========================================================
    // 打印部分数值结果
    // ==========================================================
    printf("Sample values:\n");
    printf("%-8s %-12s %-12s %-12s %-12s\n", "Cell", "x", "rho", "u", "p");
    printf("------------------------------------------------------------\n");
    
    int sample_cells[] = {0, NX/4, NX/2-1, NX/2, 3*NX/4, NX-1};
    for (int idx = 0; idx < 6; idx++) {
        int i = sample_cells[idx];
        double rho = solver.solution.uh[i][0][0];
        double rho_u = solver.solution.uh[i][0][1];
        double E = solver.solution.uh[i][0][2];
        double u = rho_u / rho;
        double p = (GAMMA- 1.0) * (E - 0.5 * rho * u * u);
        
        printf("%-8d %-12.6f %-12.6f %-12.6f %-12.6f\n", 
               i, solver.mesh.cell_centers[i], rho, u, p);
    }
    
    // ==========================================================
    // 总结
    // ==========================================================
    double elapsed = io_get_wall_time() - start_time;
    
    printf("\n");
    printf("========================================\n");
    printf("  Initialization Test Summary\n");
    printf("========================================\n");
    printf("  Status:     SUCCESS ✓\n");
    printf("  Wall time:  %.6f seconds\n", elapsed);
    printf("\n");
    printf("Next steps:\n");
    printf("  1. Visualize results:\n");
    printf("     python plot_Euler.py\n");
    printf("  2. Add spatial discretization\n");
    printf("  3. Add time integration\n");
    printf("\n");
    
    return 0;
}
