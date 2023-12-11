import math
import numpy as np
import matplotlib.pyplot as plt


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
    print(P)
    print(Q)
    return(x)

x = [0, 1, 2, 3, 4]
I = [-0.6666666667, 7.3333333333, 44.3333333333, 159.363]
n = len(x) - 1
h = x[1] - x[0]
A = [[0] * (n + 2) for i in range(n + 1)]
print(A)
for i in range(1, n):
    A[i][i - 1] = 1/h
    A[i][i] = 4/h
    A[i][i+1] = 1/h
    A[i][-1] = 3*(I[i-1] + I[i])
A[0][0] = 1
A[0][1] = 2
A[0][-1] = (5/h**2*I[0]+1/h**2*I[1])/2
A[n][n] = 1
A[n][n-1] = 2
A[n][-1] = (5/h**2*I[3]+1/h**2*I[2])/2
for a in A:
    print(a)
m = progonka(A)
print(m)
S = [0] * n
"""for i in range(n):
    S[i] = f[i] + np.poly1d([x[i]], True) * (1/h[i+1]*(f[i+1]-f[i]) - h[i+1]/2*m[i] - h[i+1]/6*(m[i+1]-m[i])) + m[i]/2*np.poly1d([x[i], x[i]], True) + 1/(6*h[i+1])*(m[i+1] - m[i])*np.poly1d([x[i], x[i], x[i]], True)"""
for i in range(n):
    S[i] = np.poly1d([-6/h**3*(I[i]-h*m[i]) + 3/h**2*(m[i+1] - m[i]), 6/h**2*(I[i]-h*m[i]) - 2/h*(m[i+1] - m[i]), m[i]])
print(S)
for i in range(n):
   print(f"{S[i][0]} + {S[i][0 + 1]} * (x - {x[i]}) + {S[i][0 + 2]} * (x - {x[i]})^2")

plt.title("Интерполяционный кубический сплайн")
x_interp = np.linspace(np.min(x), np.max(x), 5000)
#plt.plot(x, I, "o", label="Data points")
for i in range(len(x)-1):
    xi = np.linspace(x[i], x[i + 1])
    plt.plot(xi, I[i], "red", label="First line")

for i in range(len(S)):
     t = np.linspace(x[i], x[i + 1], 1000)
     plt.plot(t, S[i](t - x[i]), "green")
plt.plot(0, 0, "green", label="Cubic splain")

plt.legend()
plt.show()