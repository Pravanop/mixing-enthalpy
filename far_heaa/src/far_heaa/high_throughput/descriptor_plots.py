import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib

matplotlib.rcParams['font.size'] = 12
# Data
x = np.array([0, 0.25, 0.5, 0.75, 1])
n = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
#
# # Plotting
# fig, ax = plt.subplots(figsize=(8, 6))
#
# # Define the viridis colormap
# cmap = cm.get_cmap('viridis', len(x) - 1)
#
# for idx in range(1, len(x)):
#     y = (-x[idx] + n + 1) / abs((x[idx] - 2 * (n + 1)))
#     ax.plot(n, y, marker='o', linewidth=3, markeredgecolor='black', ms=10, alpha=0.8, color=cmap(idx - 1))
#     ax.fill_between(n, y, y2=0, alpha=0.3, color=cmap(idx - 1), label ='_nolegend_')  # Shading below the line
#
# plt.ylim([0.3, 0.5])
# plt.xlim((1, 16))
# plt.xlabel('n (order)', fontsize = 12)
# plt.ylabel('$H_c/H_a$', fontsize = 12)
# # Legend
# plt.legend(x[1:], title='x', fontsize = 12)
# plt.show()

# ha_hc = [-3, -1, -0.05, 0.05, 1, 3]
#
# cmap = cm.get_cmap('plasma', len(ha_hc) - 1)
# for idx, i in enumerate(ha_hc):
#     y = 1 / ((n + 1) ** 2)
#     y *= -(2*n+1) + n*i
#     print(y)
#     plt.plot(n, y, marker = 'o', linewidth= 2, markeredgecolor = 'black', color=cmap(idx - 1))
#
# plt.legend(ha_hc, title = 'Ha/Hc', frameon = False)
# plt.xlabel('Order of alloy (n)')
# plt.ylabel('$\delta H_n/H_n$')
# plt.savefig('delta_hn_hn.png')




# for ind, Hc in enumerate([-0.4, -0.05, 0.05, 0.4]):
#
#     for dividend in [-300, -100, -20, 20, 100, 300]:
#         Ha = dividend/100*Hc
#
#
#         delta_hn = x[idx]/((n*(n+1))**2)
#         delta_hn *= (x[idx] - 2*n - 2*Hc + (n + 1 -x[idx])*Ha)
#
#
#         ax[ind].plot(n, delta_hn, marker = 'o', linewidth= 2, markeredgecolor = 'black')
#
#     ax[ind].set_title(f'Hc = {Hc} ev/atom,Ha = y%*Hc', y = -0.08)
#     # ax[ind].set_yticks([-0.1, -0.5, ])
#
# ax[ind].legend([-300, -100, -20, 20, 100, 300], frameon = False, title = 'y', ncols = 1, bbox_to_anchor = (1.2, 3))
# plt.subplots_adjust(wspace=0, hspace=0, right = 0.85)
# fig.supxlabel('Order of alloy (n)')
# fig.supylabel('$\delta H_n$ (ev/atom)')
# plt.show()
#
#
# for idx in [1, 2, 3, 4]:
#
#     sn = (n + 1 - x[idx])*np.log((n+1 - x[idx])/(n*(n+1))) + x[idx]*np.log(x[idx]/(n+1)) + (n+1)*np.log(n)
#     sn *= - 1/(n+1)
#     Sn = -np.log(1/n)
#     plt.plot(n, sn/Sn, marker = 'o', linewidth= 2, markeredgecolor = 'black')
#     plt.ylabel('$\delta S_n/S_n$')
#     plt.xlabel('Order of alloy (n)')
#
# plt.legend(x[1:], frameon = False, title = 'x')
# plt.savefig('delta_sn_sn.png')


# def hn(n, x, ha, hc):
#     hn = x / ((n) * (n + 1)) ** 2
#     hn *= ((x - 2 * n - 2) * hc + (n + 1 - x) * ha)
#     return hn
#
# hc = np.linspace(-2, 2, 100)
# ha = np.linspace(-2, 2, 100)
# mesh = np.meshgrid(hc, ha)
# hc = mesh[0].flatten()
# ha = mesh[1].flatten()
# n = np.arange(2, 10)
# x =[0.25, 0.5, 0.75, 1]
# fig, ax = plt.subplots(8, 1, figsize = (10, 10))
# for idx, i in enumerate(n):
#     y = [hn(i, 1, hc[idx], ha[idx]) for idx in range(len(hc))]
#     y = np.array(y)
#     y = [1 if i <= 0 else 0 for i in y]
#     ax[idx].scatter(hc, ha, c = y , cmap = 'viridis')
#
# # plt.colorbar()
# plt.savefig('hn.png')

T = 1500
Hc = 0.4
Ha = -0.5

n = np.linspace(2, 15, 14)
print(n)
gn_gn = 1/((n*(n+1))**2) * ((-2*n-1)*Hc + n*Ha)+T*8.617e-5*np.log(n/(n+1))
# gn_gn /= Hc/n**2 +T*8.617e-5*np.log(1/n)
print(gn_gn)
plt.plot(n, gn_gn, marker = 'o', linewidth= 2, markeredgecolor = 'black')
plt.ylabel('$\delta G_n$')
plt.xlabel('Order of alloy (n)')
plt.savefig('delta_gn.png')
# plt.show()

