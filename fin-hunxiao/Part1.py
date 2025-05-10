import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['font.family'] = 'SimHei'

def solve_heat_conduction(nx, ny, dx, dy, dt, total_time, 
                                   rho, c, k, 
                                   T_init, T_boundary):
    """
    求解非稳态导热方程
    
    参数:
    nx, ny: x和y方向节点数
    dx, dy: 空间步长
    dt: 时间步长
    total_time: 总模拟时间
    rho: 密度
    c: 比热容
    k: 导热系数
    T_init: 初始温度
    T_boundary: 边界温度
    """
    # 计算热扩散率
    alpha = k / (rho * c)
    
    # 检查稳定性条件
    stability = alpha * dt * (1/dx**2 + 1/dy**2)
    print(f"稳定性参数: {stability} (应小于0.5)")
    if stability > 0.5:
        print("警告: 不满足显式方法稳定性条件，可能导致数值不稳定")
    
    # 初始化温度场
    T = np.ones((ny, nx)) * T_init
    
    # 设置边界条件（第一类边界）
    T[0, :] = T_boundary  # 下边界
    T[-1, :] = T_boundary  # 上边界
    T[:, 0] = T_boundary  # 左边界
    T[:, -1] = T_boundary  # 右边界
    
    # 保存各时间步的温度场
    T_history = [T.copy()]
    time_steps = [0]
    
    # 时间推进
    n_steps = int(total_time / dt)
    for step in range(1, n_steps + 1):
        T_old = T.copy()
        
        # 内部节点计算
        for i in range(1, ny-1):
            for j in range(1, nx-1):
                # 显式差分格式
                T[i, j] = T_old[i, j] + alpha * dt * (
                    (T_old[i+1, j] - 2*T_old[i, j] + T_old[i-1, j]) / dy**2 +
                    (T_old[i, j+1] - 2*T_old[i, j] + T_old[i, j-1]) / dx**2
                )
        
        # 更新边界
        T[0, :] = T_boundary  # 下边界
        T[-1, :] = T_boundary  # 上边界
        T[:, 0] = T_boundary  # 左边界
        T[:, -1] = T_boundary  # 右边界
        
        # 保存结果
        T_history.append(T.copy())
        time_steps.append(step * dt)
    
    return T_history, time_steps, alpha

def solve_tridiagonal(a, b, c, d):
    """
    求解三对角矩阵系统：Ax = d，其中A是三对角矩阵
    a: 下对角线元素
    b: 主对角线元素
    c: 上对角线元素
    d: 右侧向量
    """
    n = len(d)
    cp = np.zeros(n-1)
    dp = np.zeros(n)
    x = np.zeros(n)
    
    # 前向消元
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    
    for i in range(1, n-1):
        m = b[i] - a[i-1] * cp[i-1]
        cp[i] = c[i] / m
        dp[i] = (d[i] - a[i-1] * dp[i-1]) / m
    
    dp[n-1] = (d[n-1] - a[n-2] * dp[n-2]) / (b[n-1] - a[n-2] * cp[n-2])
    
    # 后向代入
    x[n-1] = dp[n-1]
    for i in range(n-2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i+1]
    
    return x

def solve_heat_conduction_ADI(nx, ny, dx, dy, dt, total_time, 
                              rho, c, k, T_init, T_boundary):
    """
    使用ADI方法求解二维非稳态热传导方程
    """
    # 计算热扩散率
    alpha = k / (rho * c)
    
    # 初始化温度场
    T = np.ones((ny, nx)) * T_init
    
    # 设置边界条件
    T[0, :] = T_boundary  # 下边界
    T[-1, :] = T_boundary  # 上边界
    T[:, 0] = T_boundary  # 左边界
    T[:, -1] = T_boundary  # 右边界
    
    # 保存结果
    T_history = [T.copy()]
    time_steps = [0]
    
    # 计算系数
    rx = alpha * dt / (2 * dx**2)
    ry = alpha * dt / (2 * dy**2)
    
    # 构建三对角矩阵系数
    # X方向的三对角矩阵系数
    ax = -rx * np.ones(nx-2)
    bx = (1 + 2*rx) * np.ones(nx-2)
    cx = -rx * np.ones(nx-2)
    
    # Y方向的三对角矩阵系数
    ay = -ry * np.ones(ny-2)
    by = (1 + 2*ry) * np.ones(ny-2)
    cy = -ry * np.ones(ny-2)
    
    # 时间推进
    n_steps = int(total_time / dt)
    for step in range(1, n_steps + 1):
        # 第一半步：x方向隐式，y方向显式
        T_half = np.copy(T)
        
        # 对每一行（y固定）求解x方向的隐式方程
        for i in range(1, ny-1):
            # 构建右侧向量
            rhs = np.zeros(nx-2)
            for j in range(1, nx-1):
                rhs[j-1] = T[i, j] + ry * (T[i+1, j] - 2*T[i, j] + T[i-1, j])
            
            # 修改边界条件影响
            rhs[0] += rx * T_boundary
            rhs[-1] += rx * T_boundary
            
            # 求解三对角矩阵
            T_half[i, 1:-1] = solve_tridiagonal(ax, bx, cx, rhs)
        
        # 第二半步：y方向隐式，x方向显式
        # 对每一列（x固定）求解y方向的隐式方程
        for j in range(1, nx-1):
            # 构建右侧向量
            rhs = np.zeros(ny-2)
            for i in range(1, ny-1):
                rhs[i-1] = T_half[i, j] + rx * (T_half[i, j+1] - 2*T_half[i, j] + T_half[i, j-1])
            
            # 修改边界条件影响
            rhs[0] += ry * T_boundary
            rhs[-1] += ry * T_boundary
            
            # 求解三对角矩阵
            T[1:-1, j] = solve_tridiagonal(ay, by, cy, rhs)
        
        # 保存结果
        T_history.append(T.copy())
        time_steps.append(step * dt)
    
    return T_history, time_steps, alpha



# 主函数
def cal(nx, ny, choose_method="explicit"):
    # 网格参数
    Lx, Ly = 0.2, 0.2  # 计算域尺寸，单位米
    dx, dy = Lx/nx, Ly/ny
    
    # 时间参数
    dt = 0.05  # 时间步长，秒
    total_time = 200  # 总模拟时间，秒
    
    # 初始条件和边界条件
    T_init = 500.0  # 初始温度，摄氏度
    T_boundary = 400.0  # 边界温度，摄氏度
    
    # 定义材料
    materials = [
        {"name": "A356铝合金", "rho": 2680, "c": 960, "k": 150},     # α = 1e-6
    ]
    
    # 存储结果
    all_results = []
    
    # 对每种材料进行计算
    time1 = time.time()
    for mat in materials:
        rho, c, k = mat["rho"], mat["c"], mat["k"]
        print(f"\n计算 {mat['name']}, ρ={rho}, c={c}, k={k}")


        if choose_method == "explicit":
            T_history, time_steps, alpha = solve_heat_conduction(
                nx, ny, dx, dy, dt, total_time, 
                rho, c, k, T_init, T_boundary
            )
        elif choose_method == "ADI":
            T_history, time_steps, alpha = solve_heat_conduction_ADI(
                nx, ny, dx, dy, dt, total_time, 
                rho, c, k, T_init, T_boundary
            )
        else:
            raise ValueError("选择的方法无效，请选择 'explicit' 或 'ADI'")

        mat["alpha"] = alpha
        mat["T_history"] = T_history
        mat["time_steps"] = time_steps
        all_results.append(mat)
        
        # print(f"{mat['name']} 的热扩散率: {alpha:.2e} m²/s")
    
    time2 = time.time()
    # 比较中心点温度变化
    plt.figure(figsize=(10, 6))
    center_i, center_j = ny // 2, nx // 2
    pos_1 = int(0.1 * ny)
    pos_2 = int(0.2 * ny)
    pos_3 = int(0.3 * ny)
    pos_4 = int(0.4 * ny)
    
    for mat in all_results:
        center_temps = [T[center_i, center_j] for T in mat["T_history"]]
        pos_1_temps = [T[pos_1, center_j] for T in mat["T_history"]]
        pos_2_temps = [T[pos_2, center_j] for T in mat["T_history"]]
        pos_3_temps = [T[pos_3, center_j] for T in mat["T_history"]]
        pos_4_temps = [T[pos_4, center_j] for T in mat["T_history"]]

        plt.plot(mat["time_steps"], pos_1_temps, label=f"Dy = {pos_1 * dy:.2f}")
        plt.plot(mat["time_steps"], pos_2_temps, label=f"Dy = {pos_2 * dy:.2f}")
        plt.plot(mat["time_steps"], pos_3_temps, label=f"Dy = {pos_3 * dy:.2f}")
        plt.plot(mat["time_steps"], pos_4_temps, label=f"Dy = {pos_4 * dy:.2f}")
        plt.plot(mat["time_steps"], center_temps, label=f"Dy = {center_i * dy:.2f}")
    
    plt.xlim(0, total_time)
    plt.title(f"节点数={ny}，迭代方法为{choose_method}，计算时间{time2-time1:.2f}s，不同位置温度随时间变化")
    plt.xlabel("时间(s)")
    plt.ylabel("温度(°C)")
    # plt.text(0.75, 0.5, f"alpha = {alpha:.2e} m^2/s", fontsize=12, transform=plt.gca().transAxes)
    plt.legend()
    plt.grid(True)
    plt.show()
    
    return all_results


# 运行比较
results_1 = cal(30, 30, 'explicit')
results_1 = cal(30, 30, 'ADI')
results_2 = cal(50, 50, 'explicit')
results_2 = cal(50, 50, 'ADI')
