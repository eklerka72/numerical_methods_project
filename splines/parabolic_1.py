import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
def parabolic_spline(x,y):
    x_interp=np.linspace(np.min(x), np.max(x), 50)
    y_parabolic=interp1d(x, y, kind="quadratic")
    plt.plot(x, y, "o", label="Data points")
    plt.plot(x_interp, y_parabolic(x_interp), "red", label="Parabolic spline")
    plt.legend()
    plt.show()

def quadratic_spline(x, y):
    n = len(x)
    h = np.diff(x)
    
    # Создание трехдиагональной матрицы
    A = np.zeros((n, n))
    A[0, 0] = 1
    A[n-1, n-1] = 1
    
    for i in range(1, n-1):
        A[i, i-1] = h[i-1]
        A[i, i] = 2 * (h[i-1] + h[i])
        A[i, i+1] = h[i]
    
    # Создание вектора свободных членов
    B = np.zeros(n)
    
    for i in range(1, n-1):
        B[i] = 3 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1])
    
    # Решение системы линейных уравнений
    c = np.linalg.solve(A, B)
    
    # Вычисление коэффициентов a, b и c
    a = np.zeros(n-1)
    b = np.zeros(n-1)
    d = np.zeros(n-1)
    
    for i in range(n-1):
        a[i] = (c[i+1] - c[i]) / (3 * h[i])
        b[i] = (y[i+1] - y[i]) / h[i] - (h[i] / 3) * (2 * c[i] + c[i+1])
        d[i] = y[i]
    
    return a, b, c[:-1], d

# Пример использования
x = np.array([1, 2, 3, 4, 5])
y = np.array([2, 4, 1, 6, 3])

a, b, c, d = quadratic_spline(x, y)

for i in range(len(a)):
    print(f"Уравнение сплайна {i+1}: {a[i]}x^2 + {b[i]}x + {c[i]} = 0")

print(parabolic_spline(x, y))