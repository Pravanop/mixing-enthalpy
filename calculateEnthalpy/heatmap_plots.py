import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from tqdm import tqdm
from matplotlib.patches import Rectangle, Arrow

def add_ele(composition, add_el, pD):
	temp_space = 11
	mol_space = 10
	mol_grid = np.round(np.linspace(0.0, 1 / (len(composition) + 1), mol_space), 2)
	temp_composition = composition + [add_el]
	misc_temp = {}
	temp_grid = np.linspace(200, 3200, temp_space)
	normalized_mol = np.round((mol_grid - np.min(mol_grid)) / (np.max(mol_grid) - np.min(mol_grid)), 2)
	for idx, mol in enumerate(tqdm(mol_grid, desc='Calculating miscibility temperature')):
		temp_list = []
		for temp in temp_grid:
			mol_1 = 1 - mol
			mol_ratio = [mol_1 / len(composition)] * len(composition) + [mol]
			temp_list.append(
				pD.find_decomp_products(composition=temp_composition, mol_ratio=mol_ratio, temperature=temp,
										lattice='BCC')[1])

		misc_temp[normalized_mol[idx]] = temp_list
	df = pd.DataFrame().from_dict(misc_temp)
	df = df.T
	df.columns = temp_grid
	df = df.T

	df = df.iloc[::-1]
	cmap = sns.cubehelix_palette(start=.5, rot=-.61, light=.98, dark=.35, hue=1, as_cmap=True)
	sns.set_theme(rc={'figure.figsize': (11.7, 8.27)})
	sns.set(font_scale=1.4, )
	ax = sns.heatmap(df, yticklabels=temp_grid[::-1].astype(int), xticklabels=normalized_mol, cmap=cmap, square=True,
					 cbar_kws={'label': '$E_{hull}$ eV/atom', "shrink": 0.8})
	plt.text(s='-'.join(composition), y=len(temp_grid) + 0.6, x=-2)
	plt.text(s='-'.join(temp_composition), y=len(temp_grid) + 0.6, x=len(mol_grid) + 0.1)
	ax.set_xlabel(f"${add_el}_x$")
	ax.set_ylabel("T (K)")
	plt.subplots_adjust(bottom=0.2, top=0.9)
	return ax
	# plt.savefig(f'heatmap_{"-".join(temp_composition)}.png')