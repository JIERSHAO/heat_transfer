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

# 测试代码，验证 TDMA 解法的正确性
def test_tdma():
    n = 10
    np.random.seed(42)
    
    # 生成对角元素
    a = np.random.rand(n-1) * 10   # 下对角线
    b = np.random.rand(n) * 20     # 主对角线
    c = np.random.rand(n-1) * 10   # 上对角线
    d = np.random.rand(n) * 100    # 右侧向量

    
    # 使用 TDMA 求解
    x_tdma = tdma_solver(a, b, c, d)
    
    # 使用 NumPy 求解作为对照
    A = np.zeros((n, n))
    # 填充主对角线
    np.fill_diagonal(A, b)
    # 填充上对角线
    np.fill_diagonal(A[:-1, 1:], c)
    # 填充下对角线
    np.fill_diagonal(A[1:, :-1], a)
    
    x_numpy = np.linalg.solve(A, d)
    
    # 比较结果
    print(f"TDMA 解: {x_tdma}")
    print(f"NumPy 解: {x_numpy}")
    
    return 0

if __name__ == "__main__":
    is_correct = test_tdma()
