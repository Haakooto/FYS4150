import numpy as np
import matplotlib.pyplot as plt

ke = 1.38935333e5

"""
This is a "quick" and dirty imlementation of more or less whole project in python
(Not last problem (10). really need to restructure quite a bit to do that one)

Before implementing I didnt quite see use of Particle class, and this particles are only stored as r, v and qom in PT.
Realised this was a bit awkward, because of sum_Eforce(), more precisely Eforce()
Also don't really need to use f(t) = x + iy. This can make it easier to find more elegant solution to problem above
Problem with this is storing history of positions. Can initiate Particles with full storing capability,
but then T and dt must be passed to them.
Might have some ideas for solving this, but will do when implementing in cpp, not now
This might also ease problem with time-dependant potential
"""


class Particle:
    def __init__(self, r0, v0, mass, charge):
        self.r = np.asarray(r0)
        self.v = np.asarray(v0)
        self.m = mass
        self.q = charge


class PenningTrap:
    def __init__(self, B0, V0, d, ppi=True):
        self.B0 = B0
        self.V0d = 2 * V0 * d ** -2  # 2V0 / d^2
        self.ppi = ppi  # particle-particle interactions

    def insert_particles(self, particles):
        self.N = len(particles)
        self.r0 = np.zeros((3, self.N))
        self.v0 = np.zeros((3, self.N))
        self.q = np.zeros(self.N)
        self.qom = np.zeros(self.N)
        self.P = particles

        for i, p in enumerate(particles):
            self.r0[:, i] = p.r
            self.v0[:, i] = p.v
            self.q[i] = p.q
            self.qom[i] = p.q / p.m

        self.w0 = self.B0 * self.qom
        self.wzsq = self.V0d * self.qom

    def Eforce(self, dx, dy, dz, q):
        """
        This is very akward. see longer comment
        """
        rmrj = np.asarray([dx, dy, dz])
        # rmrj = self.r - other.r
        norm = np.linalg.norm(rmrj)
        return q * rmrj * norm ** -3  # * other.q

    def sum_Eforce(self, x, y, z):
        Fse = np.zeros((self.N, 3))
        if not self.ppi:
            return Fse.T
        for i in range(self.N):
            for j in range(i):
                # This is bad
                f = -ke * self.Eforce(x[i] - x[j], y[i] - y[j], z[i] - z[j], self.q[j])
                Fse[i] += f * self.qom[i]
                Fse[j] -= f * self.qom[j]
        return Fse.T

    def advance(self, y):
        f, fv, z, zv = y
        Fsex, Fsey, Fsez = self.sum_Eforce(np.real(f), np.imag(f), np.real(z))

        df = fv
        dfv = -1j * self.w0 * fv + 0.5 * self.wzsq ** 2 * f - Fsex - 1j * Fsey
        dz = zv
        dzv = -self.wzsq ** 2 * z - Fsez
        return np.asarray([df, dfv, dz, dzv])

    def simulate(self, T, dt, method):
        y0 = np.zeros((4, self.N), dtype=complex)
        y0[0] = self.r0[0] + 1j * self.r0[1]
        y0[1] = self.v0[0] + 1j * self.v0[1]
        y0[2] = self.r0[2]
        y0[3] = self.v0[2]

        S = Solver(T, dt, y0, method)
        y, t = S.solve(self.advance)
        self.T = t
        self.r = np.zeros((len(t), self.N, 3), dtype=float)
        self.r[:, :, 0] = np.real(y[:, 0, :])
        self.r[:, :, 1] = np.imag(y[:, 0, :])
        self.r[:, :, 2] = np.real(y[:, 2, :])

        self.v = np.zeros((len(t), self.N, 3), dtype=float)
        self.v[:, :, 0] = np.real(y[:, 1, :])
        self.v[:, :, 1] = np.imag(y[:, 1, :])
        self.v[:, :, 2] = np.real(y[:, 3, :])


class Solver:
    def __init__(self, T, dt, y0, method="RK4"):
        nT = int(T / dt)
        self.t = np.linspace(0, T, nT)
        self.y = np.zeros((nT, *y0.shape), dtype=complex)
        self.y[0] = y0
        if method == "RK4":
            self.method = self.RK4
        elif method == "FEuler":
            self.method = self.Euler

    def solve(self, func):
        for i, t in enumerate(self.t[1:]):
            self.y[i + 1] = self.method(self.y[i], func, self.t[1])
        return self.y, self.t

    def Euler(self, func):
        pass

    def RK4(self, y, func, h):
        k1 = func(y)
        k2 = func(y + h * k1 / 2)
        k3 = func(y + h * k2 / 2)
        k4 = func(y + h * k3)
        return y + h / 6 * (k1 + k2 + k3 + k4)


def main():
    Penning = PenningTrap(9.65e1, 9.65e8, 1e4, ppi=True)
    p1 = Particle([0, 10, 0], [0, 1, 0], 1, 1)
    p2 = Particle([10, 0, 0], [0, 1, 0], 1, 1)
    Penning.insert_particles([p1, p2])
    Penning.simulate(10, 0.01, "RK4")

    plt.plot(Penning.r[:, 0, 0], Penning.r[:, 0, 1])
    plt.plot(Penning.r[:, 1, 0], Penning.r[:, 1, 1])
    plt.show()


if __name__ == "__main__":
    main()
