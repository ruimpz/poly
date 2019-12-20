import numpy as np
import scipy.integrate as integ
import scipy.optimize as opt
import matplotlib.pyplot as plt

plt.style.use("seaborn-paper")

class Poly():
    def __init__(self, mass, radius, luminosity, a, Y):
        self.M = mass
        self.R = radius
        self.L = luminosity
        self.X = (1 - Y)/(1 + a)
        self.Y = Y
        self.Z = a * self.X

    def __str__(self):
        return "Mass = {}\nRadius = {}\na = {}\nY = {}".format(self.M, self.R, self.Z/self.X, self.Y)

    def f(self, t, y, n):
        """Returns tupple.
        System of diff equations to solve."""
        y1, y2 = y
        return y2, -np.absolute(y1)**n - 2 * y2 / t

    def find_zero(self, t, y):
        """Returns float.
        Terminal event for differential equation solver"""
        return y[0]
    find_zero.terminal = True

    def get_poly(self, n, xf = 14, y0=[1, 0], first_step = 1e-10, max_step = .1):
        """Returns array.
        Solves polytrope equations for given conditions."""
        self.n = n
        y0s = [1- first_step**2/2, - first_step]
        sol = integ.solve_ivp(lambda t, y: self.f(t, y, n), (first_step, xf), y0s, max_step=max_step, \
                              vectorized=True, events=self.find_zero)
        data_sol = np.array([sol.t, sol.y[0], sol.y[1]])
        data_sol[1, -1] = 1e-18
        data0 = np.array([0, y0[0], y0[1]]).reshape(3, 1)
        self.poly = np.concatenate((data0, data_sol), axis=1)

    def get_params(self, G=6.67430e-8, R_k=8.314462618e7):
        """Returns array.
        Uses polytrope solution to get stellar parameters."""
        self.mu = 4/(5*self.X - self.Z + 3)
        self.xi_s, self.dtheta_s = self.poly[0, -1], self.poly[2, -1]
        self.rho_0 = - self.M * self.xi_s / (4*np.pi * self.R**3 * self.dtheta_s)
        self.P_0 = G * self.M**2 / (self.R**4 * 4 * np.pi * (self.n+1)*self.dtheta_s**2)
        self.T_0 = - self.mu * G * self.M / (R_k * self.R * (self.n+1) * self.xi_s * self.dtheta_s)
        self.rho = self.rho_0 * self.poly[1]**self.n
        self.P = self.P_0 * self.poly[1]**(self.n + 1)
        self.T = self.T_0 * self.poly[1]
        self.m = (self.poly[0]/self.xi_s)**2 * self.poly[2]/self.dtheta_s

    def get_emissivity(self):
        T, X, Z, rho = self.T, self.X, self.Z, self.rho
        T6 = T*1e-6
        a = 1.2e17*((1-X-Z)/(4*X))**2*np.exp(-100*T6**(-1/3))
        phi_a = 1-a+np.sqrt(a*(a+2))
        F1=(np.sqrt(a+2)-np.sqrt(a)) / (np.sqrt(a+2)+3*np.sqrt(a))
        F2=(1-F1) / (1+8.94e15*(X/(4-3*X))*T6**(-1/6)*np.exp(-102.65*T6**(-1/3)))
        F3 = 1 - F1 - F2
        e_0 = 2.38e6*X**2*rho*T6**(-2/3) * (1+0.0123*T6**(1/3)+0.0109*T6**(2/3)\
                                            +0.00095*T6) * np.exp(-33.80*T6**(-1/3)+0.27*rho**(1/2)*T6**(-3/2))
        self.e_pp = e_0 * phi_a * (.98*F1 + .96 * F2 + .721 * F3) / .98
        self.e_cno = 8.67e27 * Z * X * rho * T6**(-2/3)*(1 + .0027*T6**(1/3) - .00778*T6**(2/3) - .000149*T6) \
            * np.exp(-152.28 * T6**(-1/3))
        self.e = self.e_pp + self.e_cno

    def get_luminosity(self):
        self.l = 1 /(-self.xi_s**2 * self.dtheta_s) * \
            integ.simps(self.poly[0]**2 * self.poly[1]**self.n * self.M * self.e / self.L, x=self.poly[0])

    def get_data(self, n, xf = 14):
        self.get_poly(n, xf = xf)
        self.get_params()
        self.get_emissivity()
        self.get_luminosity()

    def get_grav_estimate(self):
        I = integ.simps(self.poly[0]**3 * self.poly[1]**self.n * self.poly[2], x=self.poly[0])
        return 3 * self.xi_s**3 * self.dtheta_s**2 / I + 5 - self.n

    def get_lum(self, n):
        self.get_data(n)
        return self.l

    def get_n(self, a=2.2, b=4.5):
        z = opt.root_scalar(lambda x: 1 - self.get_lum(x), bracket=(a, b), xtol = 1e-8)
        self.n = z.root

    def plot_poly(self):
        plt.plot(self.poly[0], self.poly[1], label="n = {}".format(self.n))
        plt.xlabel(r"$\xi$")
        plt.ylabel(r"$\theta$")
        plt.legend()

    def plot_T_vs_R(self):
        plt.plot(self.poly[0]/self.xi_s, self.T / self.T[0], label="Temperature")
        plt.xlabel(r"$\frac{R}{R_0}$")
        plt.ylabel(r"$\frac{T}{T_0}$")
        plt.legend()

    def plot_rho_vs_R(self):
        plt.plot(self.poly[0]/self.xi_s, self.rho / self.rho[0], label="Density")
        plt.xlabel(r"$\frac{R}{R_0}$")
        plt.ylabel(r"$\frac{\rho}{\rho_0}$")
        plt.legend()

    def plot_P_vs_R(self):
        plt.plot(self.poly[0]/self.xi_s, self.P / self.P[0], label="Pressure")
        plt.xlabel(r"$\frac{R}{R_0}$")
        plt.ylabel(r"$\frac{P}{P_0}$")
        plt.legend()

    def plot_e_vs_R(self):
        plt.plot(self.poly[0]/self.xi_s, self.e, label="Emissivity")
        plt.xlabel(r"$\frac{R}{R_0}$")
        plt.ylabel(r"$e$")
        plt.legend()

    def plot_ecno_vs_R(self):
        plt.plot(self.poly[0]/self.xi_s, self.e_cno, label=r"$e_{cno}$")
        plt.xlabel(r"$\frac{R}{R_0}$")
        plt.ylabel(r"$e_{cno}$")
        plt.legend()

    def plot_epp_vs_R(self):
        plt.plot(self.poly[0]/self.xi_s, self.e_pp, label=r"$e_{pp}$")
        plt.xlabel(r"$\frac{R}{R_0}$")
        plt.ylabel(r"$e_{pp}}$")
        plt.legend()

    def plot_L_vs_M(self):
        ls = np.zeros(self.poly[0].size)
        print(self.poly[0, :3])
        for i in range(2, ls.size):
            ls[i] = 1 /(-self.xi_s**2 * self.dtheta_s) * integ.simps(self.poly[0, :i]**2 * \
                    self.poly[1, :i]**self.n * self.M * self.e[:i] / self.L, x=self.poly[0, :i])
        plt.plot(self.m, ls)
        plt.xlabel(r"$m/M$")
        plt.ylabel(r"$\frac{L_r}{L}$")
