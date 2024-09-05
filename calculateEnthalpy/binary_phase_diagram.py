import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed


from calculateEnthalpy.helper_functions.grid_code import create_mol_grid, create_multinary
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram
import matplotlib
from scipy.interpolate import CubicSpline, splrep, BSpline, UnivariateSpline

matplotlib.rcParams.update({'font.size': 18})

def find_temperature(mol_ratio, composition, pD, genre, im_list):
    misc_T = pD.find_misc_temperature(composition=composition, mol_ratio=mol_ratio, lattice=genre, phase_flag=True, batch_tag=True, im=im_list)
    return misc_T if isinstance(misc_T, float) else 5000

def binary_diagram(composition, pD,genre):
	misc_T = []
	composition = sorted(composition)
	mol_grid = create_mol_grid(n=2, grid_size=30)
	misc_T_list = []
	mol_grid = pd.DataFrame(mol_grid)
	mol_grid = mol_grid.sort_values(by=[0])
	mol_grid = mol_grid.to_numpy()

	n_alloy = len(composition)
	all_combs = create_multinary(element_list=composition, no_comb=list(range(2, n_alloy + 1)))
	im_list = []
	for dimensionality, alloy_list in all_combs.items():
		if pD.im_flag:
			im_list += pD.get_intermetallic(alloy_list)

	num_threads = os.cpu_count()
	# print(num_threads)
	with ThreadPoolExecutor(max_workers=1) as executor:
		# Start the operations and mark each future with its mol_ratio
		misc_T_list = list(tqdm(
			executor.map(lambda i: find_temperature(i, composition, pD, genre, im_list), mol_grid),
			total=len(mol_grid)
		))

	# cs = UnivariateSpline(mol_grid[:, 0], misc_T_list)
	# cs.set_smoothing_factor(2)
	# x_range = np.arange(0, 1.01, 0.01)
	fig, ax = plt.subplots(figsize=(8, 8))
	ax.set_facecolor('#004488')
	# ax.plot(x_range, cs(x_range), label='Cubic Spline', color='black', linewidth=4)
	ax.plot(mol_grid[:, 0], misc_T_list, label='Cubic Spline', color='black', linewidth=4)
	T_m = [0, 0]
	T_m[1] = pD.tm.avg_T_melt([composition[0]], [1])
	T_m[0] = pD.tm.avg_T_melt([composition[1]], [1])
	ax.plot([0, 1], T_m, '--', color = 'red', linewidth=2)
	ax.fill(mol_grid[:, 0], misc_T_list, color='#DDAA33')

	# plt.scatter(x=mol_grid[:, 0], y=misc_T_list)
	ax.set_xlim(0,1)
	ax.set_xticks([0,1], sorted(composition, reverse=True))
	ax.set_ylabel('T (K)')
	ax.set_xlabel('x')
	ax.set_yticks(np.linspace(0,3000, 6))
	ax.set_ylim(0,max(T_m)+100)
	ax.set_title('-'.join(sorted(composition, reverse=True)) + " Phase Diagram")
	plt.subplots_adjust(left = 0.15, right=0.9)
	# plt.savefig(f'binary_diagram_{"-".join(sorted(composition, reverse=True))}.png', dpi=300)
	return ax, fig

# correction = True
# equi = False
#
# if correction:
# 	binary_file_path = "./new_phase_diagram/bokas_omegas_processed.json"
# else:
# 	binary_file_path = "../../data/output_data/bokasCorrected_bcc_1/all_lattices_binaries.json"
#
# end_member_path = "./new_phase_diagram/bokas_end_members_dict.json"
#
# pD = phaseDiagram(
# 	processed_binary_file_path=binary_file_path,
# 	end_member_file_path=end_member_path,
# 	grid_size=20,
# 	im_flag=True,
# 	correction=correction,
# 	equi_flag=equi)
# composition = ["Cr", "Ti"]
#
# binary_diagram(composition, pD, "min")[1].show()