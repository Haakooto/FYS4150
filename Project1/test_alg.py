import numpy as np
import matplotlib.pyplot as plt


f = lambda x: 100 * np.exp(-10 * x)
analytic = lambda x: 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)

N = 100001
x, h = np.linspace(0, 1, N, retstep=True)
U = analytic(x)

a = -np.ones(N)
c = -np.ones(N)
b = np.ones(N) * 2
g = h ** 2 * f(x)

for i in range(1, N):
    b[i] -= 1 / b[i - 1]
    g[i] += g[i - 1] / b[i - 1]

g[-1] = g[-1] / b[-1]
for i in range(N - 2, 0, -1):
    g[i] = (g[i] - c[i] * g[i + 1]) / b[i]

plt.plot(x, U, label="analytic")
plt.plot(x, g, label="numeric")
plt.legend()
plt.show()