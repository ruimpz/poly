import poly as poly
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

Ys, As = .27, .02857
Xs = (1 - Ys) / (1 + As)
Zs = As * Xs
Ms, Rs, Ls = 1.98919e33, 6.9699e10, 3.846e33
dM, dR, dL = 5e-6 * Ms, 5e-5*Rs, 5e-4*Ls

N = 1000

ns = np.zeros(N)
x0 = np.array([Ms, Rs, Ls])
dx = np.array([dM, dR, dL])

data = np.zeros((N, x0.size))

for i in range(N):
    data[i, :] = x0 + dx*np.random.normal(x0.size)
    M, R, L = data[i, :]
    ns[i] = poly.get_n(M, R, L, Xs, Zs)

n_avg = np.mean(ns)
n_sigma = np.sqrt( np.sum((ns - n_avg)**2) / (N - 1) )

print("n = {} +- {}".format(n_avg, n_sigma))


plt.hist(ns, bins=50, density=True)
plt.show()
