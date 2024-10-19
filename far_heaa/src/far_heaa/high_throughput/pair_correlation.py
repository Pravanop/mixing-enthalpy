import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.io.json_handler import JSONHandler
from far_heaa.math_operations.thermo_calculations import ThermoMaths
import matplotlib

matplotlib.rcParams.update({'font.size': 14})


system = 3
element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf']
file = JSONHandler.load_json(folder_path = './', file_name=f'{system}_add_ele_paths_total_10_wo_im_equi')
data = JSONHandler.load_json(folder_path='../database', file_name='bokas_omegas_processed')
mol_grid_size = 5
x = np.linspace(0, 1, mol_grid_size)

tm = ThermoMaths()
ele_temp = [tm.avg_T_melt(composition=i, mol_ratio = [0.25, 0.25]) for i in element_list]
alloy_temp = [tm.avg_T_melt(composition=i.split('-'), mol_ratio = [0.33, 0.33, 0.34]) for i in list(file.keys())]
# fig, ax = plt.subplots(2, 1, sharex=True, sharey=True)
# ax[0].axvline(np.mean(ele_temp), color='black', linestyle='dashed', linewidth=2, zorder= 1)
# ax[1].axvline(np.mean(alloy_temp), color='black', linestyle='dashed', linewidth=2, zorder=1)
# sns.histplot(ele_temp,
# 			 kde=True,
# 			 color='#BB5566',
# 			 edgecolor='black',
# 			 linewidth=2,
# 			 stat = 'percent',
# 			 ax=ax[0],
# 			 zorder=2)
# sns.histplot(alloy_temp,
# 			 kde=True,
# 			 color='#BB5566',
# 			 edgecolor='black',
# 			 stat = 'percent',
# 			 linewidth=2,
# 			 ax=ax[1],
# 			 zorder=2)
#
#
# ax[0].set_title('Element Melting Temperature Distribution')
# ax[1].set_title('Alloy Melting Temperature Distribution')
# ax[1].set_xlabel('Melting Temperature (K)')
# ax[0].set_ylabel('Frequency')
# ax[1].set_ylabel('Frequency')
# plt.grid(False)
# plt.savefig('../plots/alloying_plots/melting_temp_dist.png', dpi = 200)


# temp_melt_diff = []
# for key, value in file.items():
# 	for add_ele, temp_list in value.items():
# 		alloy_melting_temp = tm.avg_T_melt(composition=key.split('-'), mol_ratio=[1/system]*system)
# 		ele_melting_temp = tm.avg_T_melt(add_ele, [1])
#
# 		temp_diff = ele_melting_temp - alloy_melting_temp
# 		temp_melt_diff.append(temp_diff)
#
# sns.histplot(temp_melt_diff, bins=20, color = '#BB5566', kde = False, zorder = 1, stat='percent', alpha = 1)
# plt.axvline(0, color='black', linestyle='dashed', linewidth=2, zorder=0)
# plt.xlabel('Melting Temp diff between alloying element and alloy (K)')
# plt.ylabel('# paths')
# plt.savefig('../plots/alloying_plots/melting_temp_diff.png', dpi = 200)


fig, ax = plt.subplots(len(element_list), len(element_list), sharex=True, sharey=True
					   , figsize=(18, 18))

ax[0][0].set_xlim(-3000, 3000)
mol_frac = 1
temp_reduction1 = np.zeros((len(element_list), len(element_list)), dtype=object)
for i, ele1 in enumerate(element_list):
	temp_ele1 = []
	for j, ele2 in enumerate(element_list):
		if ele2 == ele1:
			if j == len(element_list) - 1:
				ax[j][i].set_xlabel(ele1)
			if i == 0:
				ax[j][i].set_ylabel(ele2)
			continue
			
		for key, value in file.items():
			comp = key.split('-')
			if ele1 in comp and ele2 in comp:
				continue
				
			if ele1 in comp:
				temp_array = value[ele2]
				if np.all(np.isnan(temp_array)):
					continue
				temp_array = np.array(temp_array)
				
				if not np.isnan(temp_array[0]) and not np.isnan(temp_array[mol_frac]):
					temp_ele1.append(temp_array[mol_frac]-temp_array[0])
		
		# temp_reduction1[i][j] = temp_ele1
		sns.kdeplot(temp_ele1, color = '#BB5566', ax=ax[j][i], fill=True, alpha=1, zorder=1)
		ax[j][i].set_xlabel(ele1)
		ax[j][i].set_ylabel(ele2)
		ax[j][i].axvline(0, color='black', linestyle='-.', linewidth=1.5, zorder=0)
		ax[j][i].axvline(np.mean(temp_ele1), color='#004488', linestyle='-', linewidth=2.5, zorder=0)
		ax[j][i].set_xticks([-2000, 0, 2000])
ax[0][0].set_yticks([])

fig.supylabel('Alloying Element')
fig.supxlabel('Element in Ternary Alloy')
plt.subplots_adjust(wspace=0, hspace=0, left=0.06, right=0.99, top=0.99, bottom=0.06)
plt.savefig(f'../plots/alloying_plots/pair-pair_{mol_frac}.png', dpi = 200)
		
		
		
		
