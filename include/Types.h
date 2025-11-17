#ifndef DG_TYPES_H
#define DG_TYPES_H

// 网格和多项式参数
#define NX 400
#define POLY_ORDER 2
#define DIM_PK (POLY_ORDER + 1)
#define NUM_GAUSS_POINTS 5
#define CFL_NUMBER 0.1
#define PI 3.14159265358979323846

// 高斯积分点和权重
typedef struct {
    double points[NUM_GAUSS_POINTS];
    double weights[NUM_GAUSS_POINTS];
} GaussQuadrature;

// 基函数及其导数
typedef struct {
    double phi[NUM_GAUSS_POINTS][DIM_PK];      // 基函数在积分点的值
    double phi_x[NUM_GAUSS_POINTS][DIM_PK];    // 基函数导数
    double phi_right[1][DIM_PK];               // 右边界基函数值
    double phi_left[1][DIM_PK];                // 左边界基函数值
    double mass_matrix[1][DIM_PK];             // 质量矩阵对角元
} BasisFunctions;

// 网格信息
typedef struct {
    double x_start;                  // 计算域起点
    double x_end;                    // 计算域终点
    double cell_size;                // 网格尺寸
    double half_cell_size;           // 半网格尺寸
    double cell_centers[NX];         // 单元中心坐标
    double bc_left;                  // 左边界条件类型
    double bc_right;                 // 右边界条件类型
    double end_time;                 // 终止时间（存储在这里方便访问）
} MeshInfo;

// 时间步进信息
typedef struct {
    double current_time;
    double end_time;
    double dt;
} TimeInfo;

// 解数组
typedef struct {
    double init[NX][DIM_PK][3]; 
    double uh[NX][DIM_PK][3];           // 当前解
    double uh_stage1[NX][DIM_PK][3];    // RK3第一步
    double uh_stage2[NX][DIM_PK][3];    // RK3第二步
    double rhs[NX][DIM_PK];          // 右端项
    double rhs_temp[NX][DIM_PK];     // 临时右端项
} SolutionArrays;

// 通量计算工作数组
typedef struct {
    double uh_boundary[NX+2][DIM_PK];    // 扩展边界数组
    double uh_gauss[NX][NUM_GAUSS_POINTS]; // 积分点处的解
    double flux[NX+1][2];                  // 界面通量
    double uh_right[NX+1][1];              // 右侧重构值
    double uh_left[NX+1][1];               // 左侧重构值
} FluxArrays;

// 完整的DG求解器结构
typedef struct {
    GaussQuadrature quad;
    BasisFunctions basis;
    MeshInfo mesh;
    TimeInfo time;
    SolutionArrays solution;
    FluxArrays flux_work;
    double exact_solution[NX][NUM_GAUSS_POINTS];
} DGSolver;

// 通量函数类型定义
typedef double (*FluxFunction)(double u);

#endif // DG_TYPES_H