import poly as poly
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

Ys, As = .27, .02857
A = As + .3
X = (1 - Ys) / (1 + A)
Z = A * X
Ms, Rs, Ls = 1.98919e33, 6.9699e10, 3.846e33
M, R, L = .934*Ms, .864*Rs, .47*Ls
dM, dR, dL = .006 * Ms, .005*Rs, .02*Ls

N = 300

ns = np.zeros(N)
x0 = np.array([M, R, L])
dx = np.array([dM, dR, dL])

data = np.zeros((N, x0.size))

for i in range(N):
    data[i, :] = x0 + dx*np.random.normal(x0.size)
    M, R, L = data[i, :]
    ns[i] = poly.get_n(M, R, L, X, Z, a = 2.4, b = 2.9)

n_avg = np.mean(ns)
n_sigma = np.sqrt( np.sum((ns - n_avg)**2) / (N - 1) )

print("n = {} +- {}".format(n_avg, n_sigma))

plt.title(r"$\alpha$ Cent B - Polytropic index distribution")
plt.hist(ns, bins=40, density=True, label="N = {}".format(N))
plt.legend()
plt.show()
