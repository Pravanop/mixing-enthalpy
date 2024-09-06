import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors
from tqdm import tqdm

from calculateEnthalpy.helper_functions.grid_code import create_multinary


def find_corners(umaps, composition):
	coords = umaps.iloc[:,:-2]
	composition_points = {}
	for i in coords.columns:
		ans = coords.loc[coords[i] == 1.0]
		cond = np.where(ans.to_numpy()[0]==1.0)[0][0]
		composition_points[composition[cond]] = (float(umaps.iloc[ans.index]['umap_x_coordinate']), float(umaps.iloc[ans.index]['umap_y_coordinate']))

	return composition_points


def use_umaps_Tmisc(n, pD, composition):
	umaps = pd.read_csv(f'./0_{n}space_umap.csv', index_col=0)
	x, y = umaps['umap_x_coordinate'], umaps['umap_y_coordinate']
	tm = []
	composition_points = find_corners(umaps, composition)
	n_alloy = len(composition)
	all_combs = create_multinary(element_list=composition, no_comb=list(range(2, n_alloy + 1)))
	im_list = []
	for dimensionality, alloy_list in all_combs.items():
		if pD.im_flag:
			im_list += pD.get_intermetallic(alloy_list)
	for ind, i in tqdm(umaps.iterrows(), desc = "Iterating through phase space"):
		mol_ratio = [i[str(j)] for j in range(1, n+1)]
		for lattice in ['min']:
			misc_T = pD.find_misc_temperature(composition=composition, lattice=lattice, mol_ratio=mol_ratio, batch_tag=True, im=im_list)
			if isinstance(misc_T, float):
				tm.append(misc_T)
			else:
				tm.append(misc_T[0])

	melt_x, melt_y = [],[]
	for ind, i in enumerate(tm):
		if i != -1:
			continue
		melt_x.append(x[ind])
		melt_y.append(y[ind])


	fig, ax = plt.subplots()
	ax.axis('off')
	s = ax.scatter(x, y, c=tm, cmap='coolwarm', label="_nolegend_")
	ax.scatter(melt_x, melt_y, c = 'black', label="_no_legend_")
	fig.colorbar(s, ax=ax, label="T_misc (K)")

	# cmap = plt.get_cmap('gist_rainbow')
	# number_of_segments = len(composition)
	# segment_colors = cmap(np.linspace(0, 1, number_of_segments))
	# color_list = [colors.rgb2hex(segment_colors[i]) for i in range(number_of_segments)]
	# # i = 0
	# # for key, value in composition_points.items():
	# # 	ax.scatter(value[0], value[1], s=20, c=color_list[i], label=key)
	# # 	# ax.text(value[0], value[1], key, ha='center', va='center')
	# # 	i += 1

	# ax.legend()

	return ax, fig
