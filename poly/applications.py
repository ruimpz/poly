from poly import Poly
import matplotlib.pyplot as plt
import numpy as np

plt.style.use("seaborn-paper")

def monte_carlo(N, Star, dx, a, b):
    x0 = Star.M, Star.R, Star.L
    for i in range(N):
        x = x0 + dx*np.random.normal(len(x0))
        Star.M, Star.R, Star.L = x
        Star.get_n(a = a, b=b)
        ns[i] = Star.n
    return ns

Ms, Rs, Ls, As, ys = 1.98919e33, 6.9699e10, 3.846e33, .02857, .27

fh1 = .22
fh2 = .3
a1 = 10**(fh1)*As
a2 = 10**(fh2)*As
M1, R1, L1 = 1.105*Ms, 1.225*Rs, 1.47*Ls
M2, R2, L2 = .934*Ms, .864*Rs, .47*Ls
dM1, dR1, dL1 = .007* Ms, .004*Rs, .05*Ls
dM2, dR2, dL2 = .006* Ms, .005*Rs, .02*Ls

Sun = Poly(Ms, Rs, Ls, As, ys)
Alpha_a = Poly(M1, R1, L1, a1, ys)
Alpha_b = Poly(M2, R2, L2, a2, ys)


#for Star in [Sun, Alpha_a, Alpha_b]:
#    Star.get_n(a = 2.5, b = 4.4)
#    Star.plot_L_vs_M()
#plt.gca().legend(("Sun", r"$\alpha$ Centauri A", r"$\alpha$ Centauri B"))
#plt.savefig("../figures/luminosities.pdf", dpi = 1000, transparent=True, bbox_inches="tight")
#plt.show()



N = 10000
ns = np.zeros(N)
#dx = np.array([dM, dR, dL])
dx1 = np.array([dM1, dR1, dL1])
dx2 = np.array([dM2, dR2, dL2])


ns = monte_carlo(N, Alpha_b, dx2, a = 2.4, b = 3)
n_avg = np.mean(ns)
n_std = np.std(ns)
print("n = {} +- {}".format(n_avg, 2*n_std))

plt.hist(ns, bins=51)
plt.xlabel("n")
plt.ylabel("Number of occurrences")
plt.savefig("../figures/alpha_b_10k.pdf", dpi = 1000, transparent=True, bbox_inches="tight")
plt.show()
