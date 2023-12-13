import math
import numpy as np
import matplotlib.pyplot as plt


def evaluate_function(function, point):
    return eval(str(function) + "(" + str(point) + ")")

def progonka(A):
    n = len(A) - 1

    Q = [0] * n
    P = [0] * n
    x = [0] * (n + 1)

    P[0] = A[0][1] / -A[0][0]
    Q[0] = -A[0][-1] / -A[0][0]

    for i in range(1, len(A) - 1):
        P[i] = A[i][1 + i] / (-A[i][i] - A[i][i - 1] * P[i - 1])
        Q[i] = (A[i][i - 1] * Q[i - 1] - A[i][-1]) / (-A[i][i] - A[i][i - 1] * P[i - 1])

    x[n] = (A[n][n - 1] * Q[n - 1] - A[n][-1]) / (-A[n][n] - A[n][n - 1] * P[n - 1])
    for i in range(n - 1, -1, -1):
        x[i] = round(P[i] * x[i + 1] + Q[i], 5)
    return(x)



def cubic_splain(x, f):
    n = len(x) - 1
    #f = []
    h = []
    m = [0] * (n + 1)
    for i in range(n + 1):
        # fi = evaluate_function(y, x[i])
        # f.append(round(fi, 7))
        if i != n:
            h.append(round(x[i + 1] - x[i], 7))
    h.append(h[-1])
    flag = 1
    for i in range(n - 1):
        if h[i] != h[i + 1]:
            flag = 0
            break
    if flag:
        m[0] = round((2 * f[0] - 5 * f[1] + 4 * f[2] - f[3]) / (h[0] ** 2), 7)
        m[n] = round((-f[n - 3] + 4 * f[n - 2] - 5 * f[n - 1] + 2 * f[n]) / (h[0] ** 2), 7)
        A = mas = [[0] * (n + 2) for i in range(n + 1)]
        for i in range(1, n):
            A[i][i - 1] = h[i - 1] / 6
            A[i][i] = (h[i - 1] + h[i]) / 3
            A[i][i + 1] = h[i] / 6
            A[i][-1] = (f[i + 1] - f[i]) / h[i] - (f[i] - f[i - 1]) / h[i - 1]
        A[0][0] = 1
        A[0][-1] = m[0]
        A[n][n] = 1
        A[n][-1] = m[n]
        m = progonka(A)
    S = [0] * n
    """for i in range(n):
        S[i] = f[i] + np.poly1d([x[i]], True) * (1/h[i+1]*(f[i+1]-f[i]) - h[i+1]/2*m[i] - h[i+1]/6*(m[i+1]-m[i])) + m[i]/2*np.poly1d([x[i], x[i]], True) + 1/(6*h[i+1])*(m[i+1] - m[i])*np.poly1d([x[i], x[i], x[i]], True)"""
    for i in range(n):
        S[i] = np.poly1d([1 / (6 * h[i + 1]) * (m[i + 1] - m[i]), m[i] / 2,
                          (1 / h[i + 1] * (f[i + 1] - f[i]) - h[i + 1] / 2 * m[i] - h[i + 1] / 6 * (m[i + 1] - m[i])),
                          f[i]])
    print(S)
    for i in range(n):
        print(
            f"{S[i][0]} + {S[i][0 + 1]} * (x - {x[i]}) + {S[i][0 + 2]} * (x - {x[i]})^2 + {S[i][0 + 3]} * (x - {x[i]})^3")

    plt.title("Интерполяционный кубический сплайн")
    x_interp = np.linspace(np.min(x), np.max(x), 5000)
    plt.plot(x, f, "o", label="Data points")
    plt.plot(x, f, "red", label="First line")

    for i in range(len(S)):
        t = np.linspace(x[i], x[i + 1], 1000)
        plt.plot(t, S[i](t - x[i]), "green")
    plt.plot(0, 0, "green", label="Cubic splain")

    plt.legend()
    plt.show()

f = np.array([1, math.sqrt(3)/2, 1/2, 0])
x = np.array([0, math.pi/6, math.pi/3, math.pi/2])
cubic_splain(x, f)

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

    
    for i in range(n-1):
        a[i] = (c[i+1] - c[i]) / (3 * h[i])
        b[i] = (y[i+1] - y[i]) / h[i] - (h[i] / 3) * (2 * c[i] + c[i+1])
        #d[i] = y[i]

    plt.title("Интерполяционный кубический сплайн")
    x_interp=np.linspace(np.min(x), np.max(x), 5000)
    plt.plot(x, y, "o", label="Data points")
    plt.plot(x, y, "red", label="First line")
    x_interp=np.linspace(np.min(x), np.max(x), 50)
    y_parabolic=interp1d(x, y, kind="quadratic")
    plt.plot(x, y, "o", label="Data points")
    plt.plot(x_interp, y_parabolic(x_interp), "green", label="Parabolic spline")
    plt.plot(0, 0, "green", label="Cubic splain")

    plt.legend()
    plt.show()

  
    return a, b, c[:-1]

def restoring(x, I):
    n = len(x) - 1
    h = x[1] - x[0]
    A = [[0] * (n + 2) for i in range(n + 1)]
    print(A)
    for i in range(1, n):
        A[i][i - 1] = 1 / h
        A[i][i] = 4 / h
        A[i][i + 1] = 1 / h
        A[i][-1] = 3 * (I[i - 1] + I[i])
    A[0][0] = 1
    A[0][1] = 2
    A[0][-1] = (5 / h ** 2 * I[0] + 1 / h ** 2 * I[1]) / 2
    A[n][n] = 1
    A[n][n - 1] = 2
    A[n][-1] = (5 / h ** 2 * I[3] + 1 / h ** 2 * I[2]) / 2
    for a in A:
        print(a)
    m = progonka(A)
    print(m)
    S = [0] * n
    """for i in range(n):
        S[i] = f[i] + np.poly1d([x[i]], True) * (1/h[i+1]*(f[i+1]-f[i]) - h[i+1]/2*m[i] - h[i+1]/6*(m[i+1]-m[i])) + m[i]/2*np.poly1d([x[i], x[i]], True) + 1/(6*h[i+1])*(m[i+1] - m[i])*np.poly1d([x[i], x[i], x[i]], True)"""
    for i in range(n):
        S[i] = np.poly1d([-6 / h ** 3 * (I[i] - h * m[i]) + 3 / h ** 2 * (m[i + 1] - m[i]),
                          6 / h ** 2 * (I[i] - h * m[i]) - 2 / h * (m[i + 1] - m[i]), m[i]])
    print(S)
    for i in range(n):
        print(f"{S[i][0]} + {S[i][0 + 1]} * (x - {x[i]}) + {S[i][0 + 2]} * (x - {x[i]})^2")

    plt.title("Интерполяционный восстанавливающий сплайн")
    x_interp = np.linspace(np.min(x), np.max(x), 5000)
    # plt.plot(x, I, "o", label="Data points")
    for i in range(len(I)):
        xi = np.linspace(x[i], x[i + 1])
        plt.plot(x[i], I[i], "o", color="pink")
        plt.plot(x[i + 1], I[i], "o", color="pink")
        plt.plot([x[i], x[i + 1]], [I[i], I[i]], "red")
    plt.plot(0, 0, "o", label="Data  points", color="pink")
    plt.plot(0, 0, "red", label="First line")

    for i in range(len(S)):
        t = np.linspace(x[i], x[i + 1], 1000)
        plt.plot(t, S[i](t - x[i]), "green")
    plt.plot(0, 0, "green", label="Restoring spline")

    plt.legend()
    plt.show()

def interpolated(x,y):
    n = len(x) - 1
    h = []
    for i in range(n):
        h.append(round(x[i + 1] - x[i], 7))
    h.append(h[-1])
    A = [[0] * (n + 1) for i in range(n + 1)]

    for i in range(1, n):
        A[i][i - 1] = (h[i + 1] / h[i]) ** 2
        A[i][i] = 1
        A[i][-1] = (h[i + 1]) ** 3 / 3 * (f[i - 1] / h[i] + 2 * (1 / h[i] + 1 / h[i + 1]) * f[i] + f[i + 1] / h[i + 1])
    A[0][0] = 1
    A[0][-1] = h[i] ** 3 / (6 * (h[i] + h[i + 1])) * (
                (2 * h[i] + 3 * h[i + 1]) / h[i] * f[i - 1] + (h[i] + h[i + 1]) * (h[i] + 3 * h[i + 1]) / (
                    h[i] ** 2 * h[i + 1]) * f[i] - f[i + 1] / h[i + 1])
    A[n][n] = 1
    A[n][-1] = h[i] ** 3 / (6 * (h[i] + h[i + 1])) * (
                (3 * h[i] + 2 * h[i + 1]) / h[i + 1] * f[i + 1] + (h[i] + h[i + 1]) * (3 * h[i] + h[i + 1]) / (
                    h[i] * h[i + 1] ** 2) * f[i] - f[i - 1] / h[i] ** 2)

    I = progonka(A)
    I = I[0:-1]

    restoring(x, I)



# Пример использования
#x = np.array([1, 2, 3, 4, 5])
#y = np.array([2, 4, 1, 6, 3])

#a, b, c = quadratic_spline(x, y)


#for i in range(len(a)):
    #print(f"Уравнение сплайна {i+1}: {a[i]}x^2 + {b[i]}x + {c[i]} = 0")

