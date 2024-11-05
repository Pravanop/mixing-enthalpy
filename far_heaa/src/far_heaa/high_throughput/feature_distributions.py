from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.math_operations.thermo_calculations import ThermoMaths
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from far_heaa.io.json_handler import JSONHandler

tm = ThermoMaths()
h_alloy_list = []
h_comp_list = []
ratio_list = []
system = 4
data = JSONHandler.load_json(folder_path='../database', file_name='bokas_omegas_processed')

element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf']

comps = list(MultinaryCombinations.create_multinary(element_list, [system]).values())[0]

for i in comps:
    composition = i.split('-')
    alloys = [j for j in element_list if j not in composition]
    for add_ele in alloys:
        h_comp = 0
        for j in list(MultinaryCombinations.create_multinary(composition, no_comb=[2]).values())[0]:
            h_comp += data[j]['BCC']
        h_comp_list.append(h_comp)
        h_alloy = 0
        for j in composition:
            h_alloy += data['-'.join(sorted([j, add_ele]))]['BCC']
        h_alloy_list.append(h_alloy)

        ratio_list.append(h_comp/h_alloy)

fig, ax = plt.subplots(1, 3)

sns.histplot(h_comp_list, ax = ax[0], stat = 'percent', color= '#BB5566', alpha = 1)
sns.histplot(h_alloy_list, ax = ax[1], stat = 'percent', color= '#004488', alpha = 1)
sns.histplot(ratio_list, ax = ax[2], stat = 'percent', color ='#DDAA33', alpha = 1)
ax[0].set_xlabel('$H_C$ (eV/atom)')
ax[1].set_xlabel('$H_A$ (eV/atom)')
ax[2].set_xlabel('$H_C/H_A$')
ax[2].set_yticks([])
ax[1].set_yticks([])
ax[2].set_ylabel('')
ax[1].set_ylabel('')
plt.subplots_adjust(hspace = 0, wspace=0)
inset_ax = inset_axes(ax[2], width="50%", height="60%", loc='upper center')
inset_ax.axvline(x = 0.5, c = 'black', linewidth = 1, linestyle = '-', alpha = 0.5, zorder= 0)
inset_ax.axvline(x = -0.5, c = 'black', linewidth = 1, linestyle = '-', alpha = 0.5, zorder = 0)
inset_ax.axvspan(-10, -0.5, color='green', alpha=0.3)
inset_ax.axvspan(0.5, 10, color='green', alpha=0.3)
sns.histplot(ratio_list, ax=inset_ax, stat='percent', color ='#DDAA33', zorder=1, alpha = 1)


inset_ax.set_xlim(-5, 5)
inset_ax.set_xticks([-5, 0, 5])
inset_ax.set_yticks([0, 8, 16])
inset_ax.set_ylabel('')

# Optionally, you can hide the y-axis of the inset if it's not necessary
# inset_ax.yaxis.set_visible(False)
plt.savefig(f'../plots/high_throughput/feature_distributions_{system}.png')







