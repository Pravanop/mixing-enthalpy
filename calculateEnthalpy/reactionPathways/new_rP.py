import itertools

import adjustText
import numpy as np
from matplotlib import pyplot as plt

from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram

def add_ele(composition, add_el, pD):
	mol_space = 10
	mol_grid = np.round(np.linspace(0.0, 1 / (len(composition) + 1), mol_space), 2)
	temp_composition = composition + [add_el]
	misc_temp = {}
	normalized_mol = np.round((mol_grid - np.min(mol_grid)) / (np.max(mol_grid) - np.min(mol_grid)), 2)

	for idx, mol in enumerate(mol_grid):
		mol_1 = 1 - mol
		mol_ratio = [mol_1 / len(composition)] * len(composition) + [mol]
		t_misc = pD.find_misc_temperature(composition=temp_composition, mol_ratio=mol_ratio, lattice='BCC')
		if isinstance(t_misc, float):
			misc_temp[normalized_mol[idx]] = t_misc
		else:
			misc_temp[normalized_mol[idx]] = 5000


	return misc_temp

def new_rP(composition, pD):
	all_pathways = list(itertools.permutations(composition, len(composition)))

	path_dict = {}
	for path in all_pathways:
		path = list(path)
		count = 0
		print(path)
		temp_path = []
		while count + 1 < len(path):

			temp = add_ele(path[:count+1], path[count+1], pD)
			if count == 0:
				temp[0.0] = pD.tm.meltT[path[:count+1][0]]
			temp_path.append(temp)
			count += 1

		path_dict['-'.join(path)] = temp_path

	texts = []
	fig, ax = plt.subplots()
	for key, value in path_dict.items():
		for idx, i in enumerate(value):

			x = np.array(list(i.keys())) + idx
			y = np.array(list(i.values()))
			ax.scatter(x[0], y[0], s=20, c='black', zorder=1)
			ax.axvline(x = x[0],color='#BBBBBB', alpha=0.3, linestyle='--', zorder=0)
			temp = key.split('-')[:idx+1]
			temp_text = ax.text(s='-'.join(sorted(temp)), x=x[0], y=y[0], zorder=2, fontsize=12)
			if temp_text not in texts:
				texts.append(temp_text)

			if idx == len(value) - 1:
				ax.scatter(x[-1], y[-1], s=20, c='black', zorder=1)
				# plt.axvline(x=x[-1], color='black')
				temp_text = ax.text(s='-'.join(sorted(composition)), x=x[-1], y=y[-1], zorder=2, fontsize=12)
				if temp_text not in texts:
					texts.append(temp_text)
			ax.plot(x,y, zorder=0)
	ax.set_ylim(-100,4000)
	ax.set_xticks([])
	ax.set_ylabel('Temperature (K)', fontsize=12)
	ax.set_xlabel('Reaction Coordinate', fontsize=12)
	adjustText.adjust_text(texts,
						   lim=500,
							force_text=0.25,
						   )
	return ax,fig
# plt.savefig('reaction_pathways.png', dpi=300)






