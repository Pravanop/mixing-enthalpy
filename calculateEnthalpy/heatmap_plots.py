import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from tqdm import tqdm
from matplotlib.patches import Rectangle, Arrow

from calculateEnthalpy.helper_functions.grid_code import create_multinary


def calculate_temp_for_mol(mol, temp_grid, composition, genre, im_list, pD, temp_composition):
	temp_list = []
	for temp in temp_grid:
		mol_1 = 1 - mol
		mol_ratio = [mol_1 / len(composition)] * len(composition) + [mol]
		
		e_above_hull = pD.find_decomp_products(
			composition=temp_composition,
			mol_ratio=mol_ratio,
			temperature=temp,
			lattice=genre,
			batch_tag=True,
			im=im_list
		)[1]
		if np.isclose(e_above_hull, 0, atol=0.001):
			temp_list.append(0)
		else:
			temp_list.append(e_above_hull)

	return temp_list


def add_ele(composition, add_el, pD, genre):
	temp_space = 11
	mol_space = 10
	mol_grid = np.round(np.linspace(0.0, 1 / (len(composition) + 1), mol_space), 2)
	temp_composition = composition + [add_el]

	temp_grid = np.linspace(200, 3200, temp_space)
	normalized_mol = np.round((mol_grid - np.min(mol_grid)) / (np.max(mol_grid) - np.min(mol_grid)), 2)
	n_alloy = len(temp_composition)
	all_combs = create_multinary(element_list=temp_composition, no_comb=list(range(2, n_alloy + 1)))
	im_list = []
	for dimensionality, alloy_list in all_combs.items():
		if pD.im_flag:
			im_list += pD.get_intermetallic(alloy_list)

	misc_temp = {}
	for idx, mol in enumerate(tqdm(mol_grid)):
		misc_temp[normalized_mol[idx]] = calculate_temp_for_mol(mol, temp_grid, composition, genre, im_list, pD,
																temp_composition)

	df = pd.DataFrame().from_dict(misc_temp)
	df = df.T
	df.columns = temp_grid
	df = df.T
	df = df.apply(lambda x: x * 1000)
	positions = []
	for col in df.columns:
		try:
			first_zero_index = df[df[col] == 0].index[0]
		except IndexError:
			first_zero_index = df.index[-1]

		positions.append(first_zero_index)

	# print(positions)
	# df = df.iloc[::-1]
	cmap = sns.cubehelix_palette(start=.5, rot=-.61, light=.98, dark=.35, hue=1, as_cmap=True)
	sns.set_theme(rc={'figure.figsize': (11.7, 8.27)})
	sns.set(font_scale=1.4, )
	fig, ax = plt.subplots()
	g = sns.heatmap(df, yticklabels=temp_grid.astype(int), cmap=cmap, square=True,
					cbar_kws={'label': '$E_{hull}$ meV/atom', "shrink": 0.8}, ax=ax)
	g.set_yticklabels(g.get_yticklabels(), rotation=0)
	ax.axes.invert_yaxis()
	ax.axhline(y=0, color='k', linewidth=3)
	ax.axhline(y=df.shape[1] + 1, color='k', linewidth=3)
	ax.axvline(x=0, color='k', linewidth=3)
	ax.axvline(x=10, color='k', linewidth=3)
	# for i, idx in enumerate(positions):
	# 	ax.plot(i + 0.5, idx + 0.5, 'bo')
	count_prev = np.where(temp_grid == positions[0])[0][0] + 1
	for idx, i in enumerate(positions):
		idx2 = np.where(temp_grid == i)[0][0] + 1
		if idx2 != count_prev:
			ax.plot([idx, idx], [count_prev, idx2], color='black', linestyle='--')
			count_prev = idx2
		ax.plot([idx, idx + 1], [idx2, idx2], color='black', linestyle='--')
	ax.text(s='-'.join(composition), y=-1, x=-2)
	ax.text(s='-'.join(temp_composition), y=-1, x=len(mol_grid) + 0.1)
	f_str = '1-x'
	# ax.set_xlabel(f"$(" + ''.join(composition) + ")_{1-x}" + f"{add_el}_x$")
	ax.set_xlabel(f"$" + f"{add_el}_x$")
	ax.set_ylabel("T (K)")
	plt.subplots_adjust(bottom=0.2, top=0.9)
	return fig
# plt.savefig(f'heatmap_{"-".join(temp_composition)}.png')
