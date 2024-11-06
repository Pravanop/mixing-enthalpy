import numpy as np
import matplotlib.pyplot as plt


x = np.linspace(0.01, 0.999, 80)


omega1, omega2, omega3, omega4 = 0.42, 0.42, 0.42, 0.42

hmix = omega1*x*(1-x) + omega2*x*(1-x) + omega3*(1-x)*x + omega4*(1-x)*(x)
entropy = -x*np.log(x) - (1-x)*np.log(1-x)

gibbs = hmix - 5500*entropy*10e-5

plt.plot(x, gibbs, linewidth=2)
plt.axhline(0, color='black', linestyle='--')
plt.xlim([0,1])
plt.savefig('./gibbs_free_energy.png')