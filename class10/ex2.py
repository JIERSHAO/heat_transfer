import numpy as np

def tdma_solver(a, b, c, d):
    """
    使用 TDMA (Thomas 算法) 求解三对角矩阵方程 Ax = d，其中 A 是三对角矩阵
    
    参数:
    a: 下对角线元素的数组 (长度为 n-1)
    b: 主对角线元素的数组 (长度为 n)
    c: 上对角线元素的数组 (长度为 n-1)
    d: 等式右侧的数组 (长度为 n)
    
    返回:
    x: 方程组的解 (长度为 n)
    """
    n = len(b)
    
    # 检查输入数组长度
    if len(a) != n-1 or len(c) != n-1 or len(d) != n:
        raise ValueError("输入数组维度不匹配")
    
    # 创建临时数组
    c_prime = np.zeros(n-1)
    d_prime = np.zeros(n)
    x = np.zeros(n)
    
    # 前向消元，转化为二对角矩阵
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]
    
    for i in range(1, n-1):
        denominator = b[i] - a[i-1] * c_prime[i-1]
        c_prime[i] = c[i] / denominator
        d_prime[i] = (d[i] - a[i-1] * d_prime[i-1]) / denominator
    
    d_prime[n-1] = (d[n-1] - a[n-2] * d_prime[n-2]) / (b[n-1] - a[n-2] * c_prime[n-2])
    
    # 回代求解
    x[n-1] = d_prime[n-1]
    for i in range(n-2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i+1]
    
    return x

def solve_convection_diffusion_tdma(nx, ny, max_iter=1000, tol=1e-6):
    """
    使用 TDMA 求解二维对流扩散方程
    基于已推导的离散方程: 3.6φ(i,j) - 1.5φ(i-1,j) - 0.5φ(i+1,j) - 1.1φ(i,j-1) - 0.5φ(i,j+1) = 0
    
    参数:
    nx, ny: x和y方向节点数
    max_iter: 最大迭代次数
    tol: 收敛容差
    """
    # 初始化phi场
    phi = np.zeros((ny, nx))
    
    # 设置边界条件
    phi[0, :] = 100  # 下边界
    phi[-1, :] = 200  # 上边界
    phi[:, 0] = 50   # 左边界
    phi[:, -1] = 300  # 右边界
    
    # 初始化内部点
    for i in range(1, ny-1):
        for j in range(1, nx-1):
            phi[i, j] = 100  # 初始猜测值
    
    # 预分配系数数组 - 基于离散方程的系数
    a_x = np.ones(nx-3) * (-1.5)  # 下对角系数 (x方向)
    b_x = np.ones(nx-2) * 3.6     # 主对角系数
    c_x = np.ones(nx-3) * (-0.5)  # 上对角系数 (x方向)
    
    a_y = np.ones(ny-3) * (-1.1)  # 下对角系数 (y方向)
    b_y = np.ones(ny-2) * 3.6     # 主对角系数
    c_y = np.ones(ny-3) * (-0.5)  # 上对角系数 (y方向)
    
    iteration = 0
    error = 1.0
    
    while iteration < max_iter and error > tol:
        phi_old = phi.copy()
        
        # 按行求解 (x方向)
        for i in range(1, ny-1):
            # 为每一行构建右侧向量
            d_x = np.zeros(nx-2)
            for j in range(1, nx-1):
                # 考虑y方向邻点影响
                d_x[j-1] = 1.1 * phi[i-1, j] + 0.5 * phi[i+1, j]
                
                # 边界条件贡献
                if j == 1:
                    d_x[j-1] += 1.5 * phi[i, 0]  # 左边界贡献
                if j == nx-2:
                    d_x[j-1] += 0.5 * phi[i, nx-1]  # 右边界贡献
            
            # 使用TDMA求解此行
            if nx > 3:  # 确保有内部点需要求解
                phi_row = tdma_solver(a_x, b_x, c_x, d_x)
                # 更新解
                for j in range(1, nx-1):
                    phi[i, j] = phi_row[j-1]
        
        # 按列求解 (y方向)
        for j in range(1, nx-1):
            # 为每一列构建右侧向量
            d_y = np.zeros(ny-2)
            for i in range(1, ny-1):
                # 考虑x方向邻点影响
                d_y[i-1] = 1.5 * phi[i, j-1] + 0.5 * phi[i, j+1]
                
                # 边界条件贡献
                if i == 1:
                    d_y[i-1] += 1.1 * phi[0, j]  # 下边界贡献
                if i == ny-2:
                    d_y[i-1] += 0.5 * phi[ny-1, j]  # 上边界贡献
            
            # 使用TDMA求解此列
            if ny > 3:  # 确保有内部点需要求解
                phi_col = tdma_solver(a_y, b_y, c_y, d_y)
                # 更新解
                for i in range(1, ny-1):
                    phi[i, j] = phi_col[i-1]
        
        # 计算误差
        error = np.max(np.abs(phi - phi_old))
        iteration += 1
        
        if iteration % 10 == 0:
            print(f"迭代 {iteration}, 误差: {error:.6e}")
    
    print(f"完成! 迭代次数: {iteration}, 最终误差: {error:.6e}")
    return phi


def solve_convection_diffusion_scipy(nx, ny):
    """使用SciPy稀疏矩阵求解"""
    from scipy import sparse
    from scipy.sparse import linalg as splinalg
    
    # 初始化phi场 
    phi = np.zeros((ny, nx))
    
    # 设置边界条件
    phi[0, :] = 100  # 下边界
    phi[-1, :] = 200  # 上边界
    phi[:, 0] = 50   # 左边界
    phi[:, -1] = 300  # 右边界
    
    # 内部节点数
    n_internal = (nx-2) * (ny-2)
    
    # 构建系数矩阵
    rows = []
    cols = []
    data = []
    
    # 右侧向量
    b = np.zeros(n_internal)
    
    # 对每个内部节点设置方程
    for i in range(1, ny-1):
        for j in range(1, nx-1):
            idx = (i-1) * (nx-2) + (j-1)  # 当前节点索引
            
            # 中心点系数
            rows.append(idx)
            cols.append(idx)
            data.append(3.6)
            
            # 左邻点
            if j > 1:
                rows.append(idx)
                cols.append(idx-1)
                data.append(-1.5)
            else:  # 左边界贡献
                b[idx] += 1.5 * phi[i, 0]
            
            # 右邻点
            if j < nx-2:
                rows.append(idx)
                cols.append(idx+1)
                data.append(-0.5)
            else:  # 右边界贡献
                b[idx] += 0.5 * phi[i, nx-1]
            
            # 下邻点
            if i > 1:
                rows.append(idx)
                cols.append(idx-(nx-2))
                data.append(-1.1)
            else:  # 下边界贡献
                b[idx] += 1.1 * phi[0, j]
            
            # 上邻点
            if i < ny-2:
                rows.append(idx)
                cols.append(idx+(nx-2))
                data.append(-0.5)
            else:  # 上边界贡献
                b[idx] += 0.5 * phi[ny-1, j]
    
    # 创建稀疏矩阵
    A = sparse.csr_matrix((data, (rows, cols)), shape=(n_internal, n_internal))
    
    # 求解系统
    phi_internal = splinalg.spsolve(A, b)
    
    # 将解重构为二维数组
    for i in range(1, ny-1):
        for j in range(1, nx-1):
            idx = (i-1) * (nx-2) + (j-1)
            phi[i, j] = phi_internal[idx]
    
    return phi


# 主程序
if __name__ == "__main__":
    # 求解问题 (使用5×5网格，包括边界点)
    phi = solve_convection_diffusion_tdma(nx=4, ny=4)
    print("最终的phi场:")
    print(np.round(phi, 2))  # 保留两位小数显示

    phi_scipy = solve_convection_diffusion_scipy(nx=4, ny=4)
    print("使用SciPy求解的phi场:")
    print(np.round(phi_scipy, 2))  # 保留两位小数显示