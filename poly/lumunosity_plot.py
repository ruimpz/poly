from poly import Poly
import matplotlib.pyplot as plt
import numpy as np

Ms, Rs, Ls, zs, ys = 1.98919e33, 6.9699e10, 3.846e33, .02857, .27

#dM, dR, dL = 5e-6 * Ms, 5e-5*Rs, 5e-4*Ls

Sun = Poly(Ms, Rs, Ls, zs, ys)
Sun.get_n()
Sun.plot_L_vs_M()
plt.show()
