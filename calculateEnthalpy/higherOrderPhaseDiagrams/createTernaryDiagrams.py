import os

import numpy as np
import matplotlib.pyplot as plt
import mpltern
from tqdm import tqdm
from calculateEnthalpy.helper_functions.grid_code import create_mol_grid
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram
from pymatgen.analysis.phase_diagram import PDPlotter
from matplotlib import tri
from scipy.interpolate import CubicSpline
import matplotlib.tri as mtri
from scipy.ndimage import gaussian_filter
"""
Example script to make a ternary diagram. Left bare for higher customizability. 
"""

correction = False
equi = False

if correction:
	binary_file_path = "../new_phase_diagram/bokas_omegas_processed.json"
else:
	binary_file_path = "../../data/output_data/all_lattices_binaries.json"

end_member_path = "../new_phase_diagram/bokas_end_members_dict.json"

pD = phaseDiagram(
	processed_binary_file_path=binary_file_path,
	end_member_file_path=end_member_path,
	grid_size=10,
	im_flag=True,
	correction=correction,
	equi_flag=equi)

composition = ['Cr', 'V', 'Ti']

temp_grid = np.linspace(200, 3000, 10).astype(int)
phase_diag_dict = pD.make_PD_comp_temp(composition=composition, temp_grid=temp_grid)

temps = list(phase_diag_dict.keys())
grid_size = 30
N_i = 100
skip = 1
sigma = 2
levels = 5

mol_grid = create_mol_grid(3, grid_size)
mol_grid = np.array(mol_grid)
mol_grid_i = create_mol_grid(3, N_i)

if equi:
	equi_str = "equi"
else:
	equi_str = "off-equi"
ele_list = '-'.join(composition)
folder_path = "../../plots/phase_diagrams"
plot_path = f"{folder_path}/{ele_list}-{equi_str}"
if not os.path.exists(plot_path):
	os.mkdir(plot_path)
t, l, r = mol_grid[:, 0], mol_grid[:, 1], mol_grid[:, 2]
# plt.show()
t_i, l_i, r_i = mol_grid_i[:, 0], mol_grid_i[:, 1], mol_grid_i[:, 2]
print(len(t),len(l),len(r))
triangulation = tri.Triangulation(t,l)
triangulation_i = tri.Triangulation(t_i,l_i)
for temp in tqdm(temps, desc="Creating phase diagrams"):
	stables = []
	# PDPlotter(phase_diag_dict[temp]).show()
	for idx, mol in enumerate(mol_grid):
		enthalpy, entropy, mol_ratio = pD.find_enthalpy_entropy_composition(composition=composition,
																   mol_ratio=mol,
																	lattice = 'BCC')

		is_stable = pD.check_stability(mol_ratio=mol_ratio,
									   temp=temp,
									   conv_hull=phase_diag_dict[temp],
									   entropy=entropy,
									   mix_enthalpy=enthalpy)
		if is_stable is not None:
			if np.isclose(is_stable[1], 0.0, atol=1e-2):
				stable = 0
			else:
				stable = 1
		else:
			stable = 0
		stables.append(stable)

	fig = plt.figure(figsize=(3,3))
	# ax = fig.add_subplot(1, 1, 1)
	ax = plt.subplot(projection="ternary")
	interp_cubic_geom = mtri.CubicTriInterpolator(triangulation, stables, kind='geom')
	stables_i = interp_cubic_geom(t_i, l_i)
	# mol_grid = np.array(mol_grid)
	# interp = RGI((t,l),stables)
	# interp_cubic_geom = tri.CubicTriInterpolator(triangulation, stables, kind='min_E')
	# stables_i = interp_cubic_geom(xi,yi)
	# ax.contourf(xi,yi,stables_i,levels=2)
	print('check')
	print(stables_i[4960])
	tricon = ax.tricontour(t, l, r, stables, levels=levels, alpha=0)
	# tricon = ax.tricontourf(t_i,l_i,r_i,stables_i,levels=levels)
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

	# path = tricon.collections[0].get_paths()[0].vertices
	# print(len(get_spiral()))
	path_i = int((levels + 1) / 2)
	try:
		#print(contour_data[path_i])
		x = contour_data[path_i][0]
		y = contour_data[path_i][1]
		x_g = gaussian_filter(x, sigma=sigma)
		y_g = gaussian_filter(y, sigma=sigma)
		# mask = np.arange(0,len(x),skip)
		# x = x[mask]
		# y = y[mask]
		p = np.linspace(0, 1, len(x))
		cs_x = CubicSpline(p, x)
		cs_y = CubicSpline(p, y)
		p_i = np.linspace(0, 1, N_i)
		x_i = cs_x(p_i)
		y_i = cs_y(p_i)
		z_i = 1 - x_i - y_i
		# ax.plot(x,y,1-x-y)
		ax.plot(x_g, y_g, 1 - x_g - y_g)
	except:
		pass
	# ax.scatter(t, l, r, c=stables)
	ax.grid()
	ax.set_tlabel(f"{composition[0]}")
	ax.set_llabel(f"{composition[1]}")
	ax.set_rlabel(f"{composition[2]}")
	plt.tight_layout()
	plt.savefig(f"{plot_path}/{temp}.png")
	# if equi == "equi":
	# 	plt.savefig(f"{folder_path}/{ele_list}/{temp}_equi.png")
	# else:
	# 	plt.savefig(f"{folder_path}/{ele_list}/{temp}.png")
