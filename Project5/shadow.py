import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time


class SparseMat:
    def __init__(self, diag, r):
        self.diag = diag
        self.r = r
        self.m = int(np.sqrt(len(diag)))

    def __getitem__(self, idx):
        if isinstance(idx, int):
            return self.diag[idx]
        else:
            i, j = idx
            if i < 0 or j < 0 or i >= self.m**2 or j >= self.m**2:
                return 0
            if i == j:
                return self.diag[i]
            elif abs(i - j) == 1:
                t = max(i, j)
                return (np.mod(t, self.m) != 0) * self.r
            else:
                return (abs(i - j) == self.m) * self.r

    def __mul__(self, x):
        b = np.zeros_like(x, complex)
        inds = [-self.m, -1, 0, 1, self.m]
        for i in range(self.m ** 2):
            for j in inds:
                b[i] += self[i, i + j] * overindex(x, i + j)
        return b

    def dense(self):
        n = self.m ** 2
        out = np.zeros((n, n), complex)
        for i in range(n):
            for j in range(n):
                out[i, j] = self[i, j]
        return out

def overindex(a, i):
    if i < 0 or i >= len(a):
        return 0
    return a[i]

def make_V(M, h, v0):
    m = M - 2
    return np.zeros((m, m))

def gauss_wave_packet(M, h, xc, yc, px, py, sx, sy):
    x = [i * h for i in range(M)]
    X, Y = np.meshgrid(x, x)

    t1x = -((X - xc) ** 2) / (2 * sx ** 2)
    t1y = -((Y - yc) ** 2) / (2 * sy ** 2)
    t2x = 1j * px * (X - xc)
    t2y = 1j * py * (Y - yc)

    u = np.exp(t1x + t1y + t2x + t2y)
    u /= np.sqrt(np.sum(np.conj(u) * u))
    return u

def make_AB(M, h, dt, V):
    m = M - 2
    a = np.ones(m * m, dtype=complex)
    b = np.ones(m * m, dtype=complex)
    r = 1j * dt / 2
    v = V.flatten()
    for k in range(m**2):
        ch = 4 * r / h ** 2 + r * v[k]
        a[k] += ch
        b[k] -= ch
    return SparseMat(a, -r), SparseMat(b, r)

def initialize(M, dt, h, v0, P):
    V = make_V(M, h, v0)
    u = gauss_wave_packet(M, h, *P)
    A, B = make_AB(M, h, dt, V)
    return A, B, u

def SOR(A, b, u=0, omega=1, tol=1e-13):
    max_iter = 1e5
    n = len(b)
    if isinstance(u, int):
        x = np.zeros_like(b)
    else:
        x = u
    max_res = 1
    iters = 0
    inds = [-A.m, -1, 1, A.m]
    while max_res > tol or iters < max_res:
        iters += 1
        print(iters)
        print(max_res)
        max_res = 0
        for i in range(n):
            sigma = 0
            for j in inds:
                sigma += A[i, i + j] * overindex(x, i + j)
            x[i] = (1 - omega) * x[i] + (omega / A[i]) * (b[i] - sigma)
            res = max(abs((A * x) - b))
            if res > max_res:
                max_res = res
    print()
    return x

def main():
    h = 0.05
    dt = 2.5e-5
    T = 0.008

    xc = 0.25
    sx = 0.005
    px = 200
    yc = 0.5
    sy = 0.05
    py = 0
    v0 = 0

    gauss_params = [xc, yc, px, py, sx, sy]

    M = int(1 / h)
    nT = int(T / dt)

    # M = 5
    # h = 1 / M

    A, B, u_ = initialize(M, dt, h, v0, gauss_params)
    u = u_[1:-1, 1:-1].flatten()

    # b = B * u
    # u = SOR(A, b, u=u, omega=1, tol=1e-13)

    probs = np.zeros((nT, (M - 2) ** 2 + 1))
    probs[:, 0] = np.linspace(0, T, nT)
    probs[0, 1:] = np.real(np.conj(u) * u)

    for i in range(1, nT):
        try:
            print(f"i = {i}")
            b = B * u
            u = SOR(A, b, u=u, omega=1, tol=1e-13)
            probs[i, 1:] = np.real(np.conj(u) * u)
        except KeyboardInterrupt:
            break

    np.save("python_solver", probs)




if __name__ == "__main__":
    main()
