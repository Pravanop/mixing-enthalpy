import json

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle, Arrow
from tqdm import tqdm
import seaborn as sns
from calculateEnthalpy.helper_functions.grid_code import create_multinary
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram
from calculateEnthalpy.helper_functions.thermo_math import thermoMaths

processed_file_path = "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/data/output_data/Ga_bcc_6/all_lattices_binaries.json"
pD = phaseDiagram(
	processed_file_path=processed_file_path,
grid_size=20)

# composition = ['Mo', 'Ta', 'Nb']
# equi = True
# if equi:
# 	mol_ratio = [1/len(composition)]*len(composition)
#
# # #Find Miscibility temperature
# print(pD.find_misc_temperature(composition=composition, mol_ratio=mol_ratio))
#
#
# #Find decomposition products at a particular temperature
# decomp_products, e_above_hull = pD.find_decomp_products(composition=composition,
# 							  mol_ratio=mol_ratio,
# 							  temperature=1100)
# print(decomp_products)
# print('E_above_hull: ', e_above_hull, " eV/atom")
#
# #Find mixing enthalpy of a multinary
# tm = thermoMaths()
# composition_dict =  dict(zip(composition, mol_ratio))
# composition_dict = {key: val for key, val in composition_dict.items() if val != 0.0}
# with open(processed_file_path, "r") as f:
# 	data = json.load(f)
# print(tm.calc_mutinary_multilattice_mix_Enthalpy(composition_dict, data))
#
# #Find mixing enthalpies of all subsets in a composition
# print(tm.find_subset_enthalpies(composition, data))

#plot heatmaps like Zhaohan with one added element
composition = ['Mo','Nb', 'Ta','Ti']
add_el = 'W'

temp_space = 10
mol_space = 10
mol_grid = np.round(np.linspace(0.0, 1/(len(composition)+1), mol_space),2)
temp_composition = composition + [add_el]
misc_temp = {}
temp_grid = np.linspace(200, 3000, temp_space)
normalized_mol = np.round((mol_grid - np.min(mol_grid))/(np.max(mol_grid) - np.min(mol_grid)),2)
for idx, mol in enumerate(tqdm(mol_grid, desc = 'Calculating miscibility temperature')):
	temp_list = []
	for temp in temp_grid:
		mol_1 = 1 - mol
		mol_ratio = [mol_1/len(composition)]*len(composition) + [mol]
		temp_list.append(pD.find_decomp_products(composition=temp_composition, mol_ratio=mol_ratio, temperature=temp)[1])

	misc_temp[normalized_mol[idx]] = temp_list
df = pd.DataFrame().from_dict(misc_temp)
df = df.T
df.columns = temp_grid
df = df.T

print(df)
# positions = [print(df[col].to_numpy()) for col in df]
positions = np.array([np.where(df[col].to_numpy() == 0)[0][0] for col in df])
positions = len(temp_grid) - positions

print(positions)
df = df.iloc[::-1]
cmap = sns.cubehelix_palette(start=.5, rot=-.61, light=.98, dark=.35, hue=1, as_cmap=True)
sns.set_theme(rc={'figure.figsize':(11.7,8.27)})
sns.set(font_scale=1.4, )
ax = sns.heatmap(df, yticklabels=temp_grid[::-1].astype(int), xticklabels=normalized_mol, cmap = cmap, square = True, cbar_kws={'label': '$E_{hull}$ eV/atom', "shrink":0.8})
plt.text(s = '-'.join(composition), y = len(temp_grid)+0.6, x = -2)
plt.text(s = '-'.join(temp_composition), y = len(temp_grid)+0.6, x = len(mol_grid)+0.1)
for idx, i in enumerate(normalized_mol):
	ax.add_patch(Rectangle(xy = (idx,positions[idx]),height=0.01,width = 1, fill=False, edgecolor='black', linestyle='--'))
ax.set_xlabel(f"${add_el}_x$")
ax.set_ylabel("T (K)")
plt.subplots_adjust(bottom = 0.2, top = 0.9)
plt.savefig('heatmap.png')
# plt.show()

