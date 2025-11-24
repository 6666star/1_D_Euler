import numpy as np
import matplotlib.pyplot as plt

def read_dg_solution(filename):
    """读取DG求解器输出的数据文件"""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            # 跳过注释行
            if line.startswith('#'):
                continue
            # 跳过空行
            if not line.strip():
                continue
            # 读取数据
            values = [float(x) for x in line.split()]
            data.append(values)
    
    data = np.array(data)
    # 返回 x, rho, u, p
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]

def read_exact_solution(filename):
    """读取精确解数据文件"""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            # 跳过注释行
            if line.startswith('#'):
                continue
            # 跳过空行
            if not line.strip():
                continue
            # 读取数据
            values = [float(x) for x in line.split()]
            data.append(values)
    
    data = np.array(data)
    # 返回 x, rho, u, p
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]

def plot_euler_results():
    """绘制欧拉方程计算结果对比图"""
    
    # 读取数据
    x_dg, rho_dg, u_dg, p_dg = read_dg_solution('output/euler_initial.dat')
    x_exact, rho_exact, u_exact, p_exact = read_exact_solution('output/euler_Exact.dat')
    
    # 创建2x2子图
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('1D Euler Equation: DG Solution vs Exact Solution', fontsize=14, fontweight='bold')
    
    # 绘制密度 rho
    axes[0, 0].plot(x_exact, rho_exact, 'k-', linewidth=2, label='Exact')
    axes[0, 0].plot(x_dg, rho_dg, 'ro', markersize=4, label='DG Solution')
    axes[0, 0].set_xlabel('x', fontsize=11)
    axes[0, 0].set_ylabel('Density (ρ)', fontsize=11)
    axes[0, 0].set_title('Density Distribution', fontsize=12)
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 绘制速度 u
    axes[0, 1].plot(x_exact, u_exact, 'k-', linewidth=2, label='Exact')
    axes[0, 1].plot(x_dg, u_dg, 'bo', markersize=4, label='DG Solution')
    axes[0, 1].set_xlabel('x', fontsize=11)
    axes[0, 1].set_ylabel('Velocity (u)', fontsize=11)
    axes[0, 1].set_title('Velocity Distribution', fontsize=12)
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # 绘制压力 p
    axes[1, 0].plot(x_exact, p_exact, 'k-', linewidth=2, label='Exact')
    axes[1, 0].plot(x_dg, p_dg, 'go', markersize=4, label='DG Solution')
    axes[1, 0].set_xlabel('x', fontsize=11)
    axes[1, 0].set_ylabel('Pressure (p)', fontsize=11)
    axes[1, 0].set_title('Pressure Distribution', fontsize=12)
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # 在第四个子图中显示误差信息
    axes[1, 1].axis('off')
    
    # 计算误差（如果网格点数相同）
    if len(x_dg) == len(x_exact):
        rho_error = np.linalg.norm(rho_dg - rho_exact) / np.linalg.norm(rho_exact)
        u_error = np.linalg.norm(u_dg - u_exact) / np.linalg.norm(u_exact)
        p_error = np.linalg.norm(p_dg - p_exact) / np.linalg.norm(p_exact)
        
        error_text = f"Relative L2 Errors:\n\n"
        error_text += f"Density (ρ):  {rho_error:.6e}\n"
        error_text += f"Velocity (u): {u_error:.6e}\n"
        error_text += f"Pressure (p): {p_error:.6e}\n"
    else:
        error_text = f"Data Information:\n\n"
        error_text += f"DG Solution points:   {len(x_dg)}\n"
        error_text += f"Exact Solution points: {len(x_exact)}\n"
    
    axes[1, 1].text(0.1, 0.5, error_text, fontsize=12, 
                    verticalalignment='center', family='monospace',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig('output/euler_comparison.png', dpi=150, bbox_inches='tight')
    print("Plot saved to output/euler_comparison.png")
    plt.show()

if __name__ == '__main__':
    try:
        plot_euler_results()
    except FileNotFoundError as e:
        print(f"Error: Cannot find data file - {e}")
        print("Please make sure the following files exist:")
        print("  - output/euler_initial.dat")
        print("  - output/euler_Exact.dat")
    except Exception as e:
        print(f"Error occurred: {e}")