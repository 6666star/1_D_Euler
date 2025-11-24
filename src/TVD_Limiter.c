//=====================TVD 限制器==================
#include "TVD_Limiter.h"
#include <math.h>
#include <string.h>

//================3X3矩阵乘法====================
 void mat_vec_mult(double mat[3][3], double vec[3], double result[3]) {
    for (int i = 0; i < 3; i++) {
        result[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            result[i] += mat[i][j] * vec[j];
        }
    }
}

//================3X3矩阵求逆====================
 void mat_inv_3x3(double mat[3][3], double inv[3][3]) {
    double det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
               - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
               + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    
    double inv_det = 1.0 / det;
    
    inv[0][0] = (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) * inv_det;
    inv[0][1] = (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) * inv_det;
    inv[0][2] = (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) * inv_det;
    inv[1][0] = (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) * inv_det;
    inv[1][1] = (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) * inv_det;
    inv[1][2] = (mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]) * inv_det;
    inv[2][0] = (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) * inv_det;
    inv[2][1] = (mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1]) * inv_det;
    inv[2][2] = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) * inv_det;
}

//========================== minmod 限制器函数===========================
 void minmod(double a[3], double b[3], double c[3], double result[3], 
                   double hx, double m) {
    for (int i = 0; i < 3; i++) {
        if (fabs(a[i]) < m * hx * hx) {
            result[i] = a[i];
        } else {
            int sign_a = (a[i] > 0) - (a[i] < 0);
            int sign_b = (b[i] > 0) - (b[i] < 0);
            int sign_c = (c[i] > 0) - (c[i] < 0);
            
            if (sign_a == sign_b && sign_a == sign_c && sign_a != 0) {
                // 取绝对值最小的
                double abs_a = fabs(a[i]);
                double abs_b = fabs(b[i]);
                double abs_c = fabs(c[i]);
                double min_abs = abs_a;
                if (abs_b < min_abs) min_abs = abs_b;
                if (abs_c < min_abs) min_abs = abs_c;
                result[i] = sign_a * min_abs;
            } else {
                result[i] = 0.0;
            }
        }
    }
}

void TVD_Limiter(DGSolver* solver, double uh[NX][DIM_PK][3]) {
    MeshInfo* mesh=&solver->mesh;
    double uhb[NX+2][DIM_PK][3];
    double uhmod[NX][DIM_PK][3];
    double deltaUR[3], deltaUL[3], deltaURM[3], deltaULM[3];
    int i, j, k;

    memset(uhmod, 0.0, sizeof(uhmod));
    memset(deltaUL, 0.0, sizeof(deltaUL));
    memset(deltaUR, 0.0, sizeof(deltaUR));
    memset(deltaULM, 0.0, sizeof(deltaULM));
    memset(deltaURM, 0.0, sizeof(deltaURM));

    for(i=0;i<NX;i++){
        for(j=0;j<3;j++){
            uhmod[i][0][j]=uh[i][0][j]; // 保留平均值不变

        }
    }

    //=================设置周期边界条件==================
    for(i = 0; i < NX; i++) {
        for(j = 0; j < DIM_PK; j++) {
            for(k = 0; k < 3; k++) {
                uhb[i+1][j][k] = uh[i][j][k];
            }
        }
    }

    if (mesh->bc_left == 1.0 && mesh->bc_right == 1.0) {
        for (j = 0; j < DIM_PK; j++) {
            for(k = 0; k < 3; k++) {
                uhb[0][j][k] = uh[NX-1][j][k];
                uhb[NX+1][j][k] = uh[0][j][k];
            }
        }
    } 
    if (mesh->bc_left == 2.0 && mesh->bc_right == 2.0) {
        for (int j = 0; j < DIM_PK; j++) {
              for(k = 0; k < 3; k++) {
                uhb[0][j][k] = uh[0][j][k];
                uhb[NX+1][j][k] = uh[NX-1][j][k];
                 }
            }
        }

   //=================================================

    for(i = 0; i < NX; i++) {
        for(j = 0; j < 3; j++) {
            deltaUR[j] = uh[i][1][j] + (2.0/3.0) * uh[i][2][j];
            deltaUL[j] = uh[i][1][j] - (2.0/3.0) * uh[i][2][j];
            deltaURM[j] = uhb[i+2][0][j] - uhb[i+1][0][j];
            deltaULM[j] = uhb[i+1][0][j] - uhb[i][0][j];
        }
    // printf("deltauR %f   %f   %f\n",deltaUR[0],deltaUR[1],deltaUR[2]);
        
        //=================特征分解==================
        double V = uh[i][0][1] / uh[i][0][0]; //速度
        double P = (GAMMA - 1) * (uh[i][0][2] - 0.5 * uh[i][0][1] * uh[i][0][1] / uh[i][0][0]); //压强 
        double C = sqrt(GAMMA * P / uh[i][0][0]); //声速
        double H = (uh[i][0][2] + P) / uh[i][0][0]; //焓
      
        double R[3][3] = {
            {1.0,      1.0,        1.0},
            {V - C,    V,          V + C},
            {H - V*C,  0.5*V*V,    H + V*C}
        };
    
        double L[3][3];
        mat_inv_3x3(R, L);

        double L_deltaUR[3], L_deltaUL[3], L_deltaURM[3], L_deltaULM[3];
        mat_vec_mult(L, deltaUR, L_deltaUR);
        mat_vec_mult(L, deltaUL, L_deltaUL);
        mat_vec_mult(L, deltaURM, L_deltaURM);
        mat_vec_mult(L, deltaULM, L_deltaULM);

        double L_deltaURM1[3], L_deltaULM1[3];
        minmod(L_deltaUR, L_deltaURM, L_deltaULM, L_deltaURM1, 
               mesh->cell_size, M);
        minmod(L_deltaUL, L_deltaURM, L_deltaULM, L_deltaULM1, 
               mesh->cell_size, M);

        double deltaURM1[3], deltaULM1[3];
        mat_vec_mult(R, L_deltaURM1, deltaURM1);
        mat_vec_mult(R, L_deltaULM1, deltaULM1);
         //printf("deltauRM1 %f   %f   %f\n",deltaULM1[0],deltaULM1[1],deltaULM1[2]);

        //======================== 更新模态系数======================
        for (int n = 0; n < 3; n++) {
            uhmod[i][1][n] = (deltaURM1[n] + deltaULM1[n]) / 2.0;
            uhmod[i][2][n] = 3.0 * (deltaURM1[n] - deltaULM1[n]) / 4.0;
        }
    }
    
    memcpy(uh, uhmod, sizeof(uhmod));


}