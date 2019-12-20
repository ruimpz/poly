import numpy as np
import poly as poly
import matplotlib.pyplot as plt


Ys, As = .27, .02857
Xs = (1 - Ys) / (1 + As)
Zs = As * Xs
Ms, Rs, Ls = 1.98919e33, 6.9699e10, 3.846e33

Ys, As = .27, .02857
A = As + .3
X = (1 - Ys) / (1 + A)
Z = A * X
Ms, Rs, Ls = 1.98919e33, 6.9699e10, 3.846e33
M, R, L = .934*Ms, .864*Rs, .47*Ls
dM, dR, dL = .006 * Ms, .005*Rs, .02*Ls


N = 100
ns = np.linspace(2, 2.8, N)
ls = ns * 0.

for i, n in enumerate(ns):
    ls[i] = poly.get_lum(n, M, R, L, X, Z)

print(ls)
plt.plot(ns, ls-1)
plt.show()
