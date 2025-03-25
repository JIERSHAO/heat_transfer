import numpy as np
from ex2_TDMA import tdma_solver

def solve_heat_conduction_tdma(T_left, T_right, T_top, n_nodes=5, max_iter=5000, tol=1e-6):
    """只使用行方向TDMA求解热传导问题"""
    # 初始化温度分布矩阵
    T = np.zeros((n_nodes, n_nodes))
    
    # 设置已知边界条件
    T[1:n_nodes-1, 0] = T_left        # 左边界
    T[1:n_nodes-1, n_nodes-1] = T_right  # 右边界
    T[0, 1:n_nodes-1] = T_top         # 顶边界
    
    # 预分配系数数组
    a = np.ones(3) * 1.0
    b = np.ones(3) * (-4.0)
    c = np.ones(3) * 1.0
    
    iterations = 0
    error = float('inf')
    
    while iterations < max_iter and error > tol:
        T_old = T.copy()
        
        # 按行求解 (x方向)
        for i in range(1, n_nodes-1):
            d = np.zeros(3)
            
            # 计算右侧向量
            for j in range(3):
                j_actual = j + 1
                d[j] = -T[i-1, j_actual] - T[i+1, j_actual]
                
                if j_actual == 1:
                    d[j] -= T[i, 0]
                if j_actual == 3:
                    d[j] -= T[i, 4]
            
            # 求解这一行
            row_solution = tdma_solver(a[:-1], b, c[:-1], d)
            
            # 更新温度
            for j in range(3):
                T[i, j+1] = row_solution[j]
        
        # 强制处理绝热底边条件
        for j in range(1, n_nodes-1):
            T[n_nodes-1, j] = T[n_nodes-2, j]
        
        # 计算误差
        error = np.max(np.abs(T - T_old))
        iterations += 1
        
    print(f"收敛需要的迭代次数: {iterations}")
    return T

# 求解并绘制结果
T = solve_heat_conduction_tdma(T_left=[37,20,5], T_right=[24,11,2], T_top=[45,35,30], n_nodes=5)
print("节点温度分布:")
print(np.round(T, 2))  # 保留一位小数显示