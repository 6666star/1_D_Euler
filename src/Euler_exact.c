#include "Euler_exact.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int solve_euler_riemann_exact(
    double rho_L, double u_L, double p_L,
    double rho_R, double u_R, double p_R,
    double gamma, double t, 
    const double *x, int n, double x0,
    double *rho, double *u, double *p)
{
    // 输入检查
    if (rho_L < 0 || rho_R < 0 || p_L < 0 || p_R < 0 || gamma <= 1.0 || t <= 0) {
        return -1;
    }
    
    // 计算音速
    double a_L = (rho_L > 0) ? sqrt(gamma * p_L / rho_L) : 0.0;
    double a_R = (rho_R > 0) ? sqrt(gamma * p_R / rho_R) : 0.0;
    
    // 求解中间压强 p_m (Newton-Raphson迭代)
    double p_m;
    if (p_L > 0 && p_R > 0) {
        // 双稀疏波初始猜测
        double exp1 = (gamma - 1.0) / (2.0 * gamma);
        double num = a_L + a_R - 0.5 * (gamma - 1.0) * (u_R - u_L);
        double den = a_L / pow(p_L, exp1) + a_R / pow(p_R, exp1);
        p_m = pow(num / den, 2.0 * gamma / (gamma - 1.0));
    } else {
        p_m = 0.5 * fmax(p_L, p_R);
    }
    
    // Newton-Raphson迭代
    if (p_m > p_L || p_m > p_R) {
        for (int iter = 0; iter < 100; iter++) {
            if (p_m <= 0) return -1;
            
            // 计算 f_L 和 f_R
            double f_L, df_L, f_R, df_R;
            double AA_L = 2.0 / ((gamma + 1.0) * rho_L);
            double BB_L = (gamma - 1.0) / (gamma + 1.0) * p_L;
            if (p_m >= p_L) {
                double dp = p_m - p_L;
                double pp = p_m + BB_L;
                double tt = sqrt(AA_L / pp);
                f_L = dp * tt;
                df_L = tt * (1.0 - 0.5 * dp / pp);
            } else {
                double pr = p_m / p_L;
                double theta = pow(pr, (gamma - 1.0) / (2.0 * gamma));
                f_L = 2.0 / (gamma - 1.0) * a_L * (theta - 1.0);
                df_L = a_L / (gamma * p_m) * theta;
            }
            
            double AA_R = 2.0 / ((gamma + 1.0) * rho_R);
            double BB_R = (gamma - 1.0) / (gamma + 1.0) * p_R;
            if (p_m >= p_R) {
                double dp = p_m - p_R;
                double pp = p_m + BB_R;
                double tt = sqrt(AA_R / pp);
                f_R = dp * tt;
                df_R = tt * (1.0 - 0.5 * dp / pp);
            } else {
                double pr = p_m / p_R;
                double theta = pow(pr, (gamma - 1.0) / (2.0 * gamma));
                f_R = 2.0 / (gamma - 1.0) * a_R * (theta - 1.0);
                df_R = a_R / (gamma * p_m) * theta;
            }
            
            double F = f_L + f_R + (u_R - u_L);
            double dF = df_L + df_R;
            double delta = F / dF;
            double p_new = p_m - delta;
            double rel_chg = fabs(delta) / fabs(0.5 * (p_m + p_new));
            p_m = p_new;
            
            if (rel_chg < 1e-8) break;
        }
    }
    
    // 计算中间速度 u_m
    double f_L, f_R;
    double AA_L = 2.0 / ((gamma + 1.0) * rho_L);
    double BB_L = (gamma - 1.0) / (gamma + 1.0) * p_L;
    if (p_m >= p_L) {
        double dp = p_m - p_L;
        double tt = sqrt(AA_L / (p_m + BB_L));
        f_L = dp * tt;
    } else {
        double theta = pow(p_m / p_L, (gamma - 1.0) / (2.0 * gamma));
        f_L = 2.0 / (gamma - 1.0) * a_L * (theta - 1.0);
    }
    
    double AA_R = 2.0 / ((gamma + 1.0) * rho_R);
    double BB_R = (gamma - 1.0) / (gamma + 1.0) * p_R;
    if (p_m >= p_R) {
        double dp = p_m - p_R;
        double tt = sqrt(AA_R / (p_m + BB_R));
        f_R = dp * tt;
    } else {
        double theta = pow(p_m / p_R, (gamma - 1.0) / (2.0 * gamma));
        f_R = 2.0 / (gamma - 1.0) * a_R * (theta - 1.0);
    }
    double u_m = 0.5 * ((u_L - f_L) + (u_R + f_R));
    
    // 计算左中间状态
    double rho_mL, a_mL, S_L, S_L_tail;
    int is_left_shock;
    if (p_m >= p_L) {
        is_left_shock = 1;
        double theta = p_m / p_L;
        double rat2 = (gamma - 1.0) / (2.0 * gamma);
        double rat3 = 1.0 - rat2;
        double rat1 = rat2 / rat3;
        rho_mL = ((rat1 + theta) / (rat1 * theta + 1.0)) * rho_L;
        a_mL = sqrt(gamma * p_m / rho_mL);
        S_L = u_L - sqrt(rat2 + rat3 * theta) * a_L;
        S_L_tail = S_L;
    } else {
        is_left_shock = 0;
        double theta = p_m / p_L;
        rho_mL = rho_L * pow(theta, 1.0 / gamma);
        a_mL = sqrt(gamma * p_m / rho_mL);
        S_L = u_L - a_L;
        S_L_tail = u_m - a_mL;
    }
    
    // 计算右中间状态
    double rho_mR, a_mR, S_R, S_R_head;
    int is_right_shock;
    if (p_m >= p_R) {
        is_right_shock = 1;
        double theta = p_m / p_R;
        double rat2 = (gamma - 1.0) / (2.0 * gamma);
        double rat3 = 1.0 - rat2;
        double rat1 = rat2 / rat3;
        rho_mR = ((rat1 + theta) / (rat1 * theta + 1.0)) * rho_R;
        a_mR = sqrt(gamma * p_m / rho_mR);
        S_R = u_m + sqrt(rat2 + rat3 / theta) * a_mR;
        S_R_head = S_R;
    } else {
        is_right_shock = 0;
        double theta = p_m / p_R;
        rho_mR = rho_R * pow(theta, 1.0 / gamma);
        a_mR = sqrt(gamma * p_m / rho_mR);
        S_R = u_R + a_R;
        S_R_head = u_m + a_mR;
    }
    
    // 在每个点计算解
    for (int i = 0; i < n; i++) {
        double xi = (x[i] - x0) / t;  // 相似变量
        
        if (xi < S_L) {
            // 左侧初始状态
            rho[i] = rho_L;
            u[i] = u_L;
            p[i] = p_L;
        } else if (is_left_shock && xi < u_m) {
            // 左激波后
            rho[i] = rho_mL;
            u[i] = u_m;
            p[i] = p_m;
        } else if (!is_left_shock && xi < S_L_tail) {
            // 左稀疏波内部
            u[i] = (2.0 / (gamma + 1.0)) * (a_L + (gamma - 1.0) / 2.0 * u_L + xi);
            double a = a_L + (gamma - 1.0) / 2.0 * (u_L - u[i]);
            rho[i] = rho_L * pow(a / a_L, 2.0 / (gamma - 1.0));
            p[i] = p_L * pow(a / a_L, 2.0 * gamma / (gamma - 1.0));
        } else if (xi < u_m) {
            // 左中间状态
            rho[i] = rho_mL;
            u[i] = u_m;
            p[i] = p_m;
        } else if (xi < S_R_head) {
            // 右中间状态
            rho[i] = rho_mR;
            u[i] = u_m;
            p[i] = p_m;
        } else if (is_right_shock && xi < S_R) {
            // 右激波前
            rho[i] = rho_mR;
            u[i] = u_m;
            p[i] = p_m;
        } else if (!is_right_shock && xi < S_R) {
            // 右稀疏波内部
            u[i] = (2.0 / (gamma + 1.0)) * (-a_R + (gamma - 1.0) / 2.0 * u_R + xi);
            double a = a_R - (gamma - 1.0) / 2.0 * (u_R - u[i]);
            rho[i] = rho_R * pow(a / a_R, 2.0 / (gamma - 1.0));
            p[i] = p_R * pow(a / a_R, 2.0 * gamma / (gamma - 1.0));
        } else {
            // 右侧初始状态
            rho[i] = rho_R;
            u[i] = u_R;
            p[i] = p_R;
        }
    }
    
    return 0;
}

void io_write_Euler_exact( double *x, double *rho, double *u, double *p,
                          int n, const char *filename)
{
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error: Cannot open file %s\n", filename);
        return;
    }

    fprintf(fp, "# 1D Euler Exact Riemann Solution\n");
    fprintf(fp, "# Columns:\n");
    fprintf(fp, "#   1: x\n");
    fprintf(fp, "#   2: rho\n");
    fprintf(fp, "#   3: u\n");
    fprintf(fp, "#   4: p\n\n");

    // 写入数值
    for (int i = 0; i < n; i++) {
        fprintf(fp, "%.10e %.10e %.10e %.10e\n",
                x[i], rho[i], u[i], p[i]);
    }

    fclose(fp);
    printf("Exact solution saved to %s\n", filename);
}
