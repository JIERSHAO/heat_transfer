import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['font.family'] = 'SimHei'

def solve_heat_absorption(nx, ny, dx, dy, dt, total_time,
                                       rho, c, k,
                                       T_init, T_inner_boundary, T_outer_boundary):
    """
    模拟砂型吸热过程
    
    参数:
    nx, ny: x和y方向节点数
    dx, dy: 空间步长
    dt: 时间步长
    total_time: 总模拟时间
    rho: 密度
    c: 比热容
    k: 导热系数
    T_init: 砂型初始温度
    T_inner_boundary: 内边界温度（铸件接触面温度）
    T_outer_boundary: 外边界温度
    """
    # 计算热扩散率和吸热系数
    alpha = k / (rho * c)
    heat_absorption_coef = np.sqrt(rho * c * k)
    
    # 检查稳定性条件
    stability = alpha * dt * (1/dx**2 + 1/dy**2)
    if stability > 0.5:
        print(f"警告: 不满足显式方法稳定性条件 (值: {stability})")
    
    # 初始化温度场
    T = np.ones((ny, nx)) * T_init
    
    # 设置内外边界条件
    # 假设中间区域是铸件，外围是砂型
    inner_start_i, inner_end_i = 2*ny//5, 3*ny//5
    inner_start_j, inner_end_j = 2*nx//4, 3*nx//5
    
    # 初始化铸件区域边界温度
    for i in range(inner_start_i, inner_end_i):
        for j in range(inner_start_j, inner_end_j):
            T[i, j] = T_inner_boundary
    
    # 设置外边界
    T[0, :] = T_outer_boundary  # 下边界
    T[-1, :] = T_outer_boundary  # 上边界
    T[:, 0] = T_outer_boundary  # 左边界
    T[:, -1] = T_outer_boundary  # 右边界
    
    # 初始化吸热量
    total_heat_absorbed = 0.0
    heat_history = [0.0]
    time_steps = [0.0]
    T_history = [T.copy()]
    
    # 时间推进
    n_steps = int(total_time / dt)
    for step in range(1, n_steps + 1):
        T_old = T.copy()
        step_heat_absorbed = 0.0
        
        # 内部砂型区域计算
        for i in range(1, ny-1):
            for j in range(1, nx-1):
                # 跳过铸件内部区域
                if (inner_start_i < i < inner_end_i-1 and 
                    inner_start_j < j < inner_end_j-1):
                    continue
                
                # 砂型区域计算
                if not (inner_start_i <= i < inner_end_i and 
                       inner_start_j <= j < inner_end_j):
                    T[i, j] = T_old[i, j] + alpha * dt * (
                        (T_old[i+1, j] - 2*T_old[i, j] + T_old[i-1, j]) / dy**2 +
                        (T_old[i, j+1] - 2*T_old[i, j] + T_old[i, j-1]) / dx**2
                    )
        
        # 计算传热量，因为铸件温度较高，始终流向砂型
        q_in = 0.0
        for i in range(inner_start_i, inner_end_i):
            for j in range(inner_start_j, inner_end_j):
                # 下边界点
                if i == inner_start_i:
                    q_in += k * (T_old[i, j] - T_old[i-1, j]) / dy
                
                # 上边界点
                if i == inner_end_i-1:
                    q_in += k * (T_old[i, j] - T_old[i+1, j]) / dy
                
                # 左边界点
                if j == inner_start_j:
                    q_in += k * (T_old[i, j] - T_old[i, j-1]) / dx
                
                # 右边界点
                if j == inner_end_j-1:
                    q_in += k * (T_old[i, j] - T_old[i, j+1]) / dx
                
                # 累加吸热量
                step_heat_absorbed += q_in * dt
                
                # 保持内边界温度
                T[i, j] = T_inner_boundary

        # 累加吸热量
        total_heat_absorbed += step_heat_absorbed    

        # 保持外边界温度
        T[0, :] = T_outer_boundary
        T[-1, :] = T_outer_boundary
        T[:, 0] = T_outer_boundary
        T[:, -1] = T_outer_boundary
        
        # 每隔一定步数保存结果
        T_history.append(T.copy())
        heat_history.append(total_heat_absorbed)
        time_steps.append(step * dt)

    return T, total_heat_absorbed, heat_history, time_steps, heat_absorption_coef, T_history

def compare_different_materials():
    # 网格参数
    nx, ny = 100, 100
    Lx, Ly = 1, 1  # 计算域尺寸，单位米
    dx, dy = Lx/(nx-1), Ly/(ny-1)
    
    # 时间参数
    dt = 10  # 时间步长，秒
    total_time = 3000  # 总模拟时间，秒
    
    # 温度条件
    T_init = 25.0  # 砂型初始温度，摄氏度
    T_inner = 400.0  # 内边界温度（铸件温度），摄氏度
    T_outer = 25.0  # 外边界温度，摄氏度
    
    # 不同的材料参数组合，变化ρ、c、λ，观察吸热系数与吸热量的关系
    materials = [
        {"name": "粘土砂", "rho": 1500, "c": 1200, "k": 1},
        {"name": "树脂砂", "rho": 1600, "c": 1000, "k": 0.7},
        {"name": "水玻璃砂", "rho": 1700, "c": 1100, "k": 1.1},
        {"name": "测试", "rho": 2000, "c": 1000, "k": 1.5},
    ]
    
    # 存储结果
    results = []
    
    # 对每种材料进行计算
    for mat in materials:
        rho, c, k = mat["rho"], mat["c"], mat["k"]
        print(f"\n计算 {mat['name']}, ρ={rho}, c={c}, k={k}")
        
        final_T, total_heat, heat_history, time_steps, heat_coef, T_history = solve_heat_absorption(
            nx, ny, dx, dy, dt, total_time,
            rho, c, k, T_init, T_inner, T_outer
        )
        
        print(f"{mat['name']} 的吸热系数: {heat_coef:.2f}, 总吸热量: {total_heat:.2f} J")
        
        # 存储结果
        mat["heat_absorption_coef"] = heat_coef
        mat["total_heat_absorbed"] = total_heat
        mat["heat_history"] = heat_history
        mat["time_steps"] = time_steps
        mat["final_temp"] = final_T
        mat["T_history"] = T_history
        results.append(mat)
    
    # 绘制吸热系数与总吸热量的关系
    plt.figure(figsize=(10, 6))
    coefs = [mat["heat_absorption_coef"] for mat in results]
    heats = [mat["total_heat_absorbed"] for mat in results]
    
    plt.scatter(coefs, heats, s=100)
    for i, mat in enumerate(results):
        plt.annotate(mat["name"], (coefs[i], heats[i]), fontsize=12)
    
    # 拟合直线
    z = np.polyfit(coefs, heats, 1)
    p = np.poly1d(z)
    plt.plot(coefs, p(coefs), "r--", linewidth=2)
    
    plt.title("砂型吸热量与吸热系数的关系")
    plt.xlabel("吸热系数 √(ρcλ)")
    plt.ylabel("吸热量 (J)")
    plt.grid(True)
    plt.show()
    
    # 绘制不同材料的吸热历史
    plt.figure(figsize=(10, 6))
    for mat in results:
        plt.plot(mat["time_steps"], mat["heat_history"], 
                 label=f"{mat['name']}")
    
    plt.title("不同砂型材料的吸热过程")
    plt.xlim(0, total_time)
    plt.ylim(0, max([max(mat["heat_history"]) for mat in results]) * 1.1)
    plt.xlabel("时间 (s)")
    plt.ylabel("吸热量 (J)")
    plt.legend()
    plt.grid(True)
    plt.show()

    # 绘制温度变化
    plt.figure(figsize=(10, 6))

    # 只显示第一种材料的结果
    mat = results[0]

    # 计算(0.39m, 0.5m)位置对应的网格索引
    pos_i = int(0.5 * ny)  # y方向
    pos_j = int(0.39 * nx)  # x方向

    # 提取该位置的温度随时间变化
    temp_history = []
    for temp_field in mat["T_history"]:
        temp_history.append(temp_field[pos_i, pos_j])

    # 绘制温度曲线
    plt.plot(mat["time_steps"], temp_history, 
            label=f"{mat['name']}在位置(0.39m, 0.5m)处")

    plt.title("位置(0.39m, 0.5m)处的温度随时间变化")
    plt.xlabel("时间 (s)")
    plt.ylabel("温度 (°C)")
    plt.legend()
    plt.grid(True)
    plt.show()

# 运行比较
results = compare_different_materials()
