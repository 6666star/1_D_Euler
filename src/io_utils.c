// ==================== io_utils.c ====================
#include "io_utils.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/time.h>
#endif

// ============================================================
// 1. 输出物理量 - 推荐用于可视化和分析
// ============================================================
void io_write_solution(const DGSolver* solver, const char* filename) {
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error: Cannot open file %s\n", filename);
        return;
    }
    
    // 写入头信息
    fprintf(fp, "# 1D Euler DG Solution - Physical Variables\n");
    fprintf(fp, "# Time: %.15e\n", solver->time.current_time);
    fprintf(fp, "# NX=%d, Polynomial Order=%d\n", NX, POLY_ORDER);
    fprintf(fp, "# gamma=%.15e\n", GAMMA);
    fprintf(fp, "#\n");
    fprintf(fp, "# Columns:\n");
    fprintf(fp, "#   1: x (cell center)\n");
    fprintf(fp, "#   2: rho (density)\n");
    fprintf(fp, "#   3: u (velocity)\n");
    fprintf(fp, "#   4: p (pressure)\n");
    fprintf(fp, "#   5: E (total energy)\n");
    fprintf(fp, "#   6: e (internal energy per unit mass)\n");
    fprintf(fp, "#   7: Ma (Mach number)\n");
    fprintf(fp, "#\n");
    
    // 输出每个单元中心的物理量
    for (int i = 0; i < NX; i++) {
        // 单元平均值 = 第0个模态系数
        double rho = solver->solution.uh[i][0][0];
        double rho_u = solver->solution.uh[i][0][1];
        double E = solver->solution.uh[i][0][2];
        
        // 计算物理量
        double u = rho_u / rho;
        double p = (GAMMA - 1.0) * (E - 0.5 * rho * u * u);
        double e =( p / (rho * GAMMA - 1.0));  // 内能
        double c = sqrt(GAMMA * p / rho);       // 声速
        double Ma = fabs(u) / c;                // 马赫数
        
        fprintf(fp, "%20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e\n",
                solver->mesh.cell_centers[i], rho, u, p, E, e, Ma);
    }
    
    fclose(fp);
    printf("Solution written to %s (t=%.6e)\n", filename, solver->time.current_time);
}

// ============================================================
// 2. 输出模态系数 - 用于重启动和精确存储
// ============================================================
void io_write_modal_coefficients(const DGSolver* solver, const char* filename) {
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error: Cannot open file %s\n", filename);
        return;
    }
    
    fprintf(fp, "# 1D Euler DG Solution - Modal Coefficients\n");
    fprintf(fp, "# Time: %.15e\n", solver->time.current_time);
    fprintf(fp, "# NX=%d, DIM_PK=%d\n", NX, DIM_PK);
    fprintf(fp, "#\n");
    fprintf(fp, "# Format: cell_id  x_center  ");
    for (int j = 0; j < DIM_PK; j++) {
        fprintf(fp, "rho_%d  ", j);
    }
    for (int j = 0; j < DIM_PK; j++) {
        fprintf(fp, "rhou_%d  ", j);
    }
    for (int j = 0; j < DIM_PK; j++) {
        fprintf(fp, "E_%d  ", j);
    }
    fprintf(fp, "\n");
    
    for (int i = 0; i < NX; i++) {
        fprintf(fp, "%5d %20.15e ", i, solver->mesh.cell_centers[i]);
        
        // rho的所有模态
        for (int j = 0; j < DIM_PK; j++) {
            fprintf(fp, "%25.16e ", solver->solution.uh[i][j][0]);
        }
        
        // rho*u的所有模态
        for (int j = 0; j < DIM_PK; j++) {
            fprintf(fp, "%25.16e ", solver->solution.uh[i][j][1]);
        }
        
        // E的所有模态
        for (int j = 0; j < DIM_PK; j++) {
            fprintf(fp, "%25.16e ", solver->solution.uh[i][j][2]);
        }
        
        fprintf(fp, "\n");
    }
    
    fclose(fp);
    printf("Modal coefficients written to %s\n", filename);
}

// ============================================================
// 3. 输出高分辨率数据 - 用于精细分析和绘图
// ============================================================
void io_write_high_resolution(const DGSolver* solver, const char* filename) {
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error: Cannot open file %s\n", filename);
        return;
    }
    
    fprintf(fp, "# 1D Euler DG Solution - High Resolution Output\n");
    fprintf(fp, "# Time: %.15e\n", solver->time.current_time);
    fprintf(fp, "# Output at all Gauss quadrature points\n");
    fprintf(fp, "#\n");
    fprintf(fp, "# Columns: x  rho  u  p  E  Ma\n");
    fprintf(fp, "#\n");
    
    for (int i = 0; i < NX; i++) {
        for (int g = 0; g < NUM_GAUSS_POINTS; g++) {
            // 计算高斯点的物理坐标
            double x = solver->mesh.cell_centers[i] + 
                       solver->mesh.half_cell_size * solver->quad.points[g];
            
            // 重构高斯点处的守恒变量
            double rho = 0.0, rho_u = 0.0, E = 0.0;
            for (int j = 0; j < DIM_PK; j++) {
                rho   += solver->solution.uh[i][j][0] * solver->basis.phi[g][j];
                rho_u += solver->solution.uh[i][j][1] * solver->basis.phi[g][j];
                E     += solver->solution.uh[i][j][2] * solver->basis.phi[g][j];
            }
            
            // 计算物理量
            double u = rho_u / rho;
            double p = (GAMMA - 1.0) * (E - 0.5 * rho * u * u);
            double c = sqrt(GAMMA * p / rho);
            double Ma = fabs(u) / c;
            
            fprintf(fp, "%20.15e %20.15e %20.15e %20.15e %20.15e %20.15e\n",
                    x, rho, u, p, E, Ma);
        }
    }
    
    fclose(fp);
    printf("High resolution solution written to %s\n", filename);
}

// ============================================================
// 4. 输出诊断信息 - 检查守恒性和数值稳定性
// ============================================================
void io_write_diagnostics(const DGSolver* solver, const char* filename, int append) {
    FILE* fp = fopen(filename, append ? "a" : "w");
    if (fp == NULL) {
        printf("Error: Cannot open file %s\n", filename);
        return;
    }
    
    if (!append) {
        fprintf(fp, "# DG Solver Diagnostics\n");
        fprintf(fp, "# Columns: time  total_mass  total_momentum  total_energy  ");
        fprintf(fp, "min_rho  max_rho  min_p  max_p  max_Ma\n");
    }
    
    double total_mass = 0.0;
    double total_momentum = 0.0;
    double total_energy = 0.0;
    double min_rho = 1e30, max_rho = -1e30;
    double min_p = 1e30, max_p = -1e30;
    double max_Ma = 0.0;
    
    for (int i = 0; i < NX; i++) {
        double rho = solver->solution.uh[i][0][0];
        double rho_u = solver->solution.uh[i][0][1];
        double E = solver->solution.uh[i][0][2];
        
        double u = rho_u / rho;
        double p = (GAMMA - 1.0) * (E - 0.5 * rho * u * u);
        double c = sqrt(GAMMA * p / rho);
        double Ma = fabs(u) / c;
        
        double dx = solver->mesh.cell_size;
        total_mass += rho * dx;
        total_momentum += rho_u * dx;
        total_energy += E * dx;
        
        if (rho < min_rho) min_rho = rho;
        if (rho > max_rho) max_rho = rho;
        if (p < min_p) min_p = p;
        if (p > max_p) max_p = p;
        if (Ma > max_Ma) max_Ma = Ma;
    }
    
    fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
            solver->time.current_time, total_mass, total_momentum, total_energy,
            min_rho, max_rho, min_p, max_p, max_Ma);
    
    fclose(fp);
}

// ============================================================
// 5. 输出单个单元的详细信息 - 调试用
// ============================================================
void io_write_cell_detail(const DGSolver* solver, int cell_id, const char* filename) {
    if (cell_id < 0 || cell_id >= NX) {
        printf("Error: Invalid cell_id %d\n", cell_id);
        return;
    }
    
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error: Cannot open file %s\n", filename);
        return;
    }
    
    fprintf(fp, "# Detailed output for cell %d\n", cell_id);
    fprintf(fp, "# Time: %.15e\n", solver->time.current_time);
    fprintf(fp, "# Cell center: %.15e\n", solver->mesh.cell_centers[cell_id]);
    fprintf(fp, "#\n");
    
    fprintf(fp, "# Modal coefficients:\n");
    for (int j = 0; j < DIM_PK; j++) {
        fprintf(fp, "# Mode %d: rho=%.15e  rho*u=%.15e  E=%.15e\n", j,
                solver->solution.uh[cell_id][j][0],
                solver->solution.uh[cell_id][j][1],
                solver->solution.uh[cell_id][j][2]);
    }
    
    fprintf(fp, "#\n");
    fprintf(fp, "# Values at quadrature points:\n");
    fprintf(fp, "# xi  x  rho  u  p  E\n");
    
    for (int g = 0; g < NUM_GAUSS_POINTS; g++) {
        double xi = solver->quad.points[g];
        double x = solver->mesh.cell_centers[cell_id] + 
                   solver->mesh.half_cell_size * xi;
        
        double rho = 0.0, rho_u = 0.0, E = 0.0;
        for (int j = 0; j < DIM_PK; j++) {
            rho   += solver->solution.uh[cell_id][j][0] * solver->basis.phi[g][j];
            rho_u += solver->solution.uh[cell_id][j][1] * solver->basis.phi[g][j];
            E     += solver->solution.uh[cell_id][j][2] * solver->basis.phi[g][j];
        }
        
        double u = rho_u / rho;
        double p = (GAMMA - 1.0) * (E - 0.5 * rho * u * u);
        
        fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e\n",
                xi, x, rho, u, p, E);
    }
    
    fclose(fp);
    printf("Cell %d details written to %s\n", cell_id, filename);
}

// ============================================================
// 6. 计时函数
// ============================================================
double io_get_wall_time(void) {
#if defined(_WIN32)
    LARGE_INTEGER freq, counter;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart / (double)freq.QuadPart;
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
#endif
}

// ============================================================
// 7. 打印运行时信息
// ============================================================
void io_print_progress(const DGSolver* solver, double wall_time) {
    double progress = solver->time.current_time / solver->mesh.end_time * 100.0;
    
    printf("\n");
    printf("========================================\n");
    printf("  Time: %.6e / %.6e (%.1f%%)\n", 
           solver->time.current_time, solver->mesh.end_time, progress);
    printf("  dt:   %.6e\n", solver->time.dt);
    printf("  Wall time: %.3f s\n", wall_time);
    printf("========================================\n");
}