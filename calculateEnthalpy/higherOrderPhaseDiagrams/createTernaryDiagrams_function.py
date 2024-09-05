import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import matplotlib.pyplot as plt
import mpltern
from tqdm import tqdm
from calculateEnthalpy.helper_functions.grid_code import create_mol_grid, create_multinary
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram
from pymatgen.analysis.phase_diagram import PDPlotter
from matplotlib import tri
from scipy.interpolate import CubicSpline
import matplotlib.tri as mtri
from scipy.ndimage import gaussian_filter
"""
Example script to make a ternary diagram. Left bare for higher customizability. 
"""


def check_mol_stability(mol, composition, genre, T, im_list, pD):
	enthalpy, entropy, mol_ratio = pD.find_enthalpy_entropy_composition(
		composition=composition,
		mol_ratio=mol,
		lattice=genre,
		temperature=T
	)

	is_stable = pD.check_stability(
		mol_ratio=mol_ratio,
		temp=T,
		conv_hull=pD.make_convex_hull(
			composition=composition,
			temperature=T,
			batch_tag=True,
			im=im_list
		),
		entropy=entropy,
		mix_enthalpy=enthalpy
	)

	if is_stable is not None:
		if np.isclose(is_stable[1], 0.0, atol=1e-2):
			return 0
		else:
			return 1
	else:
		return 0

def ternary_diagram(composition, T, pD, genre):
	grid_size = 20
	N_i = 100
	skip = 1
	sigma = 1
	levels = 5

	mol_grid = create_mol_grid(3, grid_size)
	mol_grid = np.array(mol_grid)
	mol_grid_i = create_mol_grid(3, N_i)
	t, l, r = mol_grid[:, 0], mol_grid[:, 1], mol_grid[:, 2]
	t_i, l_i, r_i = mol_grid_i[:, 0], mol_grid_i[:, 1], mol_grid_i[:, 2]
	# print(len(t),len(l),len(r))
	triangulation = tri.Triangulation(t,l)
	triangulation_i = tri.Triangulation(t_i,l_i)
	stables = []
	n_alloy = len(composition)
	all_combs = create_multinary(element_list=composition, no_comb=list(range(2, n_alloy + 1)))
	im_list = []
	for dimensionality, alloy_list in all_combs.items():
		if pD.im_flag:
			im_list += pD.get_intermetallic(alloy_list)
	num_threads = os.cpu_count()  # Example: twice the number of CPU cores for I/O-bound tasks

	# Parallel processing with ThreadPoolExecutor
	with ThreadPoolExecutor(max_workers=1) as executor:
		# Use map to parallelize the execution and ensure order is preserved
		stables = list(
			tqdm(executor.map(lambda mol: check_mol_stability(mol, composition, genre, T, im_list, pD), mol_grid),
				 total=len(mol_grid)))

	fig = plt.figure()
	ax = plt.subplot(projection="ternary")
	interp_cubic_geom = mtri.CubicTriInterpolator(triangulation, stables, kind='geom')
	stables_i = interp_cubic_geom(t_i, l_i)
	print('check')
	print(stables_i[4960])
	tricon = ax.tricontour(t, l, r, stables, levels=levels, alpha=0)
	ax.scatter(t,l,r,c=stables)
	contour_data = []
	i = 0
	for collection in tricon.collections:
		print(i)
		i += 1

		for path in collection.get_paths():
			# Extract the vertices of the contour
			v = path.vertices
			x_data = v[:, 0]
			y_data = v[:, 1]
			ternary_coords = ax.transProjection.inverted().transform(np.c_[x_data, y_data])

			contour_data.append(
				(ternary_coords[:, 0], ternary_coords[:, 1], 1 - ternary_coords[:, 0] - ternary_coords[:, 1]))

	path_i = int((levels + 1) / 2)
	try:
		x = contour_data[path_i][0]
		y = contour_data[path_i][1]
		x_g = gaussian_filter(x, sigma=sigma)
		y_g = gaussian_filter(y, sigma=sigma)
		p = np.linspace(0, 1, len(x))
		cs_x = CubicSpline(p, x)
		cs_y = CubicSpline(p, y)
		p_i = np.linspace(0, 1, N_i)
		x_i = cs_x(p_i)
		y_i = cs_y(p_i)
		z_i = 1 - x_i - y_i
		ax.plot(x_g, y_g, 1 - x_g - y_g)
	except:
		pass
	ax.grid()
	ax.set_tlabel(f"{composition[0]}")
	ax.set_llabel(f"{composition[1]}")
	ax.set_rlabel(f"{composition[2]}")
	plt.tight_layout()
	return ax, fig
