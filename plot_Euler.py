import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def check_file_exists(filename):
    if not os.path.exists(filename):
        print(f"ERROR: {filename} not found!")
        print("Please run './bin/test_init' first")
        sys.exit(1)

def plot_initial_conditions():
    print("Plotting 1D Euler initial conditions...")
    
    # 检查文件
    check_file_exists('output/euler_initial.dat')
    
    # 读取数据
    data = np.loadtxt('output/euler_initial.dat')
    
    x = data[:, 0]
    rho = data[:, 1]
    u = data[:, 2]
    p = data[:, 3]
    E = data[:, 4]
    
    # 创建图形
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('1D Euler Equations - Initial Conditions (Sod Shock Tube)', 
                 fontsize=16, fontweight='bold')
    
    # 密度
    axes[0, 0].plot(x, rho, 'b-', linewidth=2.5, label='DG P%d' % (len(data[0]) // 3 - 2))
    axes[0, 0].axvline(x=0.5, color='r', linestyle='--', linewidth=1.5, 
                       alpha=0.7, label='Discontinuity')
    axes[0, 0].axhline(y=1.0, color='gray', linestyle=':', alpha=0.5)
    axes[0, 0].axhline(y=0.125, color='gray', linestyle=':', alpha=0.5)
    axes[0, 0].set_xlabel('x', fontsize=13)
    axes[0, 0].set_ylabel('Density ρ', fontsize=13)
    axes[0, 0].set_title('Density Profile', fontsize=14)
    axes[0, 0].grid(True, alpha=0.3, linestyle='--')
    axes[0, 0].legend(fontsize=11)
    axes[0, 0].set_ylim([-0.05, 8])
    
    # 速度
    axes[0, 1].plot(x, u, 'g-', linewidth=2.5, label='DG Solution')
    axes[0, 1].axvline(x=0.5, color='r', linestyle='--', linewidth=1.5, 
                       alpha=0.7, label='Discontinuity')
    axes[0, 1].axhline(y=0.0, color='gray', linestyle=':', alpha=0.5)
    axes[0, 1].set_xlabel('x', fontsize=13)
    axes[0, 1].set_ylabel('Velocity u', fontsize=13)
    axes[0, 1].set_title('Velocity Profile', fontsize=14)
    axes[0, 1].grid(True, alpha=0.3, linestyle='--')
    axes[0, 1].legend(fontsize=11)
    axes[0, 1].set_ylim([-0.1, 30])
    
    # 压力
    axes[1, 0].plot(x, p, 'r-', linewidth=2.5, label='DG Solution')
    axes[1, 0].axvline(x=0.5, color='r', linestyle='--', linewidth=1.5, 
                       alpha=0.7, label='Discontinuity')
    axes[1, 0].axhline(y=1.0, color='gray', linestyle=':', alpha=0.5)
    axes[1, 0].axhline(y=0.1, color='gray', linestyle=':', alpha=0.5)
    axes[1, 0].set_xlabel('x', fontsize=13)
    axes[1, 0].set_ylabel('Pressure p', fontsize=13)
    axes[1, 0].set_title('Pressure Profile', fontsize=14)
    axes[1, 0].grid(True, alpha=0.3, linestyle='--')
    axes[1, 0].legend(fontsize=11)
    axes[1, 0].set_ylim([-0.05, 1100])
    
    # 总能量
    axes[1, 1].plot(x, E, 'm-', linewidth=2.5, label='DG Solution')
    axes[1, 1].axvline(x=0.5, color='r', linestyle='--', linewidth=1.5, 
                       alpha=0.7, label='Discontinuity')
    axes[1, 1].set_xlabel('x', fontsize=13)
    axes[1, 1].set_ylabel('Total Energy E', fontsize=13)
    axes[1, 1].set_title('Energy Profile', fontsize=14)
    axes[1, 1].grid(True, alpha=0.3, linestyle='--')
    axes[1, 1].legend(fontsize=11)
    
    plt.tight_layout()
    plt.savefig('output/euler_initial.png', dpi=150, bbox_inches='tight')
    print("✓ Figure saved: output/euler_initial.png")
    plt.show()

if __name__ == '__main__':
    plot_initial_conditions()