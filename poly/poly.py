import numpy as np
import scipy.integrate as integ
import scipy.optimize as opt
import matplotlib.pyplot as plt


def f(t, y, n):
    """Returns tupple.

    System of diff equations to solve."""
    y1, y2 = y
    return y2, -np.absolute(y1)**n - 2 * y2 / t


def find_zero(t, y):
    """Returns float.

    Terminal event for differential equation solver"""
    return y[0]
find_zero.terminal = True


def polytropes(xf, n, y0, first_step = 1e-10, max_step = .1):
    """Returns array.

    Solves polytrope equations for given conditions."""
    y0s = [1- first_step**2/2, - first_step]
    sol = integ.solve_ivp(lambda t, y: f(t, y, n), (first_step, xf), y0s, max_step=max_step, vectorized=True, events=find_zero)
    data_sol = np.array([sol.t, sol.y[0], sol.y[1]])
    data_sol[1, -1] = 1e-18
    data0 = np.array([0, y0[0], y0[1]]).reshape(3, 1)
    return np.concatenate((data0, data_sol), axis=1)

def get_params(data, M, R, n, X, Z, G=6.67430e-8, R_k=8.314462618e7):
    """Returns array.

    Uses polytrope solution to get stellar parameters."""
    mu = 4/(5*X - Z + 3)
    xi_s, dtheta_s = data[0, -1], data[2, -1]
    rho_0 = - M * xi_s / (4*np.pi * R**3 * dtheta_s)
    P_0 = G * M**2 / (R**4 * 4 * np.pi * (n+1)*dtheta_s**2)
    T_0 = - mu * G * M / (R_k * R * (n+1) * xi_s * dtheta_s)
    rho = rho_0 * data[1]**n
    P = P_0 * data[1]**(n + 1)
    T = T_0 * data[1]
    return rho, P, T


def get_emissivity(rho, T, X, Z):
    T6 = T*1e-6
    a = 1.2e17*((1-X-Z)/(4*X))**2*np.exp(-100*T6**(-1/3))
    phi_a = 1-a+np.sqrt(a*(a+2))
    F1=(np.sqrt(a+2)-np.sqrt(a)) / (np.sqrt(a+2)+3*np.sqrt(a))
    F2=(1-F1) / (1+8.94e15*(X/(4-3*X))*T6**(-1/6)*np.exp(-102.65*T6**(-1/3)))
    F3 = 1 - F1 - F2
    e_0 = 2.38e6*X**2*rho*T6**(-2/3) * (1+0.0123*T6**(1/3)+0.0109*T6**(2/3)\
    +0.00095*T6) * np.exp(-33.80*T6**(-1/3)+0.27*rho**(1./2.)*T6**(-3/2))
    e_pp = e_0 * phi_a * (.98*F1 + .96 * F2 + .721 * F3) / .98
    e_cno = 8.67e27 * Z * X * rho * T6**(-2/3)*(1 + .0027*T6**(1/3) - .00778*T6**(2/3) - .000149*T6) \
        * np.exp(-152.28 * T6**(-1/3))
    return e_pp + e_cno

def get_luminosity(data, M, n, e, L0):
    xi_s, dtheta_s = data[0, -1], data[2, -1]
    return 1 /(-xi_s**2 * dtheta_s) * integ.simps(data[0]**2 * data[1]**n * M * e / L0, x=data[0])



def get_lum(n, M, R, L0, X, Z, xf = 18, y0=[1, 0]):
    poli = polytropes(xf, n, y0)
    rho, P, T = get_params(poli, M, R, n, X, Z)
    e = get_emissivity(rho, T, X, Z)
    return get_luminosity(poli, M, n, e, L0)


def get_n(M, R, L, X, Z, a=3, b=3.5):
    z = opt.root_scalar(lambda x: 1 - get_lum(x, M, R, L, X, Z), bracket=(a, b), xtol = 1e-8)
    return z.root

#poli = politropos(16, z.root, [1, 0])
#eta_s, dtheta_s =poli[0, -1],poli[2, -1]
#g = 1 / eta_s**3 / dtheta_s**2 * integ.simps(-poli[0]**3 *poli[1]**z.root *poli[2], x=poli[0])


#n = 3.19
#y0 = [1, 0]
#xf = 16
#poli = politropos(xf, n, y0)
#rho, P, T = get_params(poli, M, R, n, X, Z)
#e = get_emissivity(rho, T, X, Z)
#L = get_luminosity(poli, M, n, e, L0)

#plt.plot(poli[0] / poli[0, -1], e, label="e")
#plt.plot(poli[0] / poli[0, -1], P,label=r"$P$")
#plt.legend()
#plt.show()
#plt.plot(poli[0] / poli[0, -1], T, label =r"$T$")
#plt.legend()
#plt.show()


#for n in np.arange(0, 5.5, .5):
#    plt.plot(poli[0], poli[1], label="n = {}".format(n))
#    print(poli[0, -1])
#plt.legend()
#plt.show()
