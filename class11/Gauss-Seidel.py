import numpy as np
import matplotlib.pyplot as plt


def Gauss_Seidel_FDMA(coff, b, x0, error_bound, max_iter):
    """
    Gauss-Seidel迭代法求解五对角矩阵线性方程组Ax=b

    :param coff: 系数矩阵A
    :param b: 常数向量b
    :param x0: 初始值
    :param error_bound: 误差范围
    :param max_iter: 最大迭代次数
    :return: 解向量x
    """
    x = x0.copy()
    for _ in range(max_iter):
        x_old = x.copy()
        for i in range(len(b)):
            sum1 = np.dot(coff[i, :i], x[:i])
            sum2 = np.dot(coff[i, i+1:], x_old[i+1:])
            x[i] = (b[i] - sum1 - sum2) / coff[i, i]

        # 计算误差
        error_vector = np.abs(x - x_old)
        error_sum = np.sum(error_vector)
        if error_sum < error_bound:
            break

    return x

coff = np.zeros((50, 50))
a = 12*np.ones(50)
b = -2*np.ones(49)
c = 1*np.ones(48)

# 使用NumPy的diag函数填充
coff = np.zeros((50, 50))
# np.fill_diagonal(coff, a)                   # 主对角线
# np.fill_diagonal(coff[:-1, 1:], b)          # 主对角线上方第一条
# np.fill_diagonal(coff[1:, :-1], b)          # 主对角线下方第一条
# np.fill_diagonal(coff[:-2, 2:], c)          # 主对角线上方第二条
# np.fill_diagonal(coff[2:, :-2], c)          # 主对角线下方第二条

x = np.zeros(50)  # 初始值
error_bound = 1e-10  # 误差范围
# 创建一个空列表来保存每次迭代的误差
error_history = []

for num in range(1000):
    x_old = x.copy()  # 备份当前值
    for i in range(50):
        if i == 0:
            x[i] = 1/12*(5 + 2*x_old[i+1] - x_old[i+2])
        elif i == 1:
            x[i] = 1/12*(5 + 2*x[i-1] + 2*x_old[i+1] - x_old[i+2])
        elif i == 49:
            x[i] = 1/12*(5 + 2*x[i-1] - x[i-2])
        elif i == 48:
            x[i] = 1/12*(5 + 2*x[i-1] + 2*x[i+1] - x[i-2])
        else:
            x[i] = 1/12*(5 + 2*x[i-1] + 2*x_old[i+1] - x[i-2] - x_old[i+2])  # 迭代公式


    # 计算误差
    error_vector = np.abs(x - x_old)
    error_sum = np.sum(error_vector)
    error_history.append(error_sum)
    if error_sum < error_bound:
        print(f"迭代次数: {num+1}")
        break

# 输出所有迭代的误差值
print("\n每次迭代的误差值:")
for idx, err in enumerate(error_history):
    print(f"迭代 {idx+1}: {err:.8e}")

# 可视化误差变化 (可选)
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(error_history)+1), error_history)
plt.semilogy(range(1, len(error_history)+1), error_history)
plt.xlabel('number of iterations')
plt.ylabel('error')
plt.grid(True)
plt.show()

print("最终结果:", x)
