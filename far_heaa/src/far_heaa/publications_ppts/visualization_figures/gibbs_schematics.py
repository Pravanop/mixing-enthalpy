import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.size'] = 14

D = [0.0052, 0.0096, 0.0132, 0.0159, 0.0177, 0.0185, 0.0184, 0.0174, 0.0155, 0.0128, 0.0094, 0.0056, 0.0016, -0.0015, 0.0]
ABC = [0.0052, 0.0036, 0.002, 0.0003, -0.0015, -0.0033, -0.0051, -0.0068, -0.0086, -0.0103, -0.0118, -0.0132, -0.0142, -0.0146, -0.013]

AB = [0.0052, 0.0054, 0.0057, 0.0063, 0.0072, 0.0083, 0.0098, 0.0116, 0.0138, 0.0164, 0.0196, 0.0235, 0.0283, 0.0347, 0.0453]
CD = [0.0052, 0.0054, 0.0057, 0.0063, 0.0072, 0.0083, 0.0098, 0.0116, 0.0138, 0.0164, 0.0196, 0.0235, 0.0283, 0.0347, 0.0453]

AB_CD = AB[::-1][:-1] + CD
trans = [0.0103, 0.007, 0.0064, 0.0068, 0.0079, 0.0095, 0.0114, 0.0138, 0.0164, 0.0195, 0.0229, 0.0268, 0.0314, 0.037, 0.0453]
ABC_D = ABC[::-1][:-1] + D
x = np.linspace(0, 1, len(AB_CD))
# x = np.linspace(0, 1, len(trans))
x1 = np.linspace(0, 1, len(ABC))/3
x2 = np.linspace(0.333, 1, len(D))
fig, ax = plt.subplots(2, 1)
# ax[0].plot(x, ABC_D, linewidth = 2, color = '#225555', alpha = 0.8, marker = 'o', markeredgecolor = 'black')
ax[0].plot(x1[:], ABC[::-1][:], linewidth = 2, color = '#225555', alpha = 0.8, marker = 'o', markeredgecolor = 'black')
ax[0].plot(x2, D, linewidth = 2, color = '#225555', alpha = 0.8, marker = 'o', markeredgecolor = 'black')
ax[1].plot(x, AB_CD, linewidth = 2, color = '#225555', alpha = 0.8, marker = 'o', markeredgecolor = 'black')
# ax.plot(x, trans, linewidth = 4, color = '#225555', alpha = 0.9, marker = 'o', markeredgecolor = 'black', ms =12)
# ax.set_ylabel('$G_{mix}$ (eV/atom)')
fig.supylabel('$G_{mix}$ (eV/atom)')
ax[0].set_xticks(ticks=[0, 0.33, 1], labels = ['ACD', 'ABCD', 'B'])
ax[1].set_xticks(ticks=[0, 0.5, 1], labels = ['AD', 'ABCD', 'BC'])
# ax.set_xticks(ticks=[0, 1], labels = ['DAC', 'CAB'])
plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.1)
plt.tight_layout()
plt.savefig('/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/far_heaa/src/far_heaa/publications_ppts/visualization_figures/gibbs_schematic_ABCD.png', dpi= 200)