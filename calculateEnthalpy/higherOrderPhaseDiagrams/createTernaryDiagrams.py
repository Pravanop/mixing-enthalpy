import os

import numpy as np
import matplotlib.pyplot as plt
import mpltern
from tqdm import tqdm
from calculateEnthalpy.helper_functions.grid_code import create_mol_grid
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram
from pymatgen.analysis.phase_diagram import PDPlotter

"""
Example script to make a ternary diagram. Left bare for higher customizability. 
"""

correction = False
equi = True

if correction:
	binary_file_path = "../new_phase_diagram/bokas_omegas_processed.json"
else:
	binary_file_path = "../../data/output_data/bokasCorrected_bcc_1/all_lattices_binaries.json"

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
mol_grid = create_mol_grid(3, grid_size)
mol_grid = np.array(mol_grid)
ele_list = '-'.join(composition)
folder_path = "../../plots/phase_diagrams"
if not os.path.exists(f"{folder_path}/{ele_list}"):
	os.mkdir(f"{folder_path}/{ele_list}")
t, l, r = mol_grid[:, 0], mol_grid[:, 1], mol_grid[:, 2]
# plt.show()
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

	fig, ax = plt.subplots(1, 1)
	ax = plt.subplot(projection="ternary")
	mol_grid = np.array(mol_grid)
	ax.scatter(t, l, r,c= stables)
	ax.grid()
	ax.set_tlabel(f"{ele_list[0]}")
	ax.set_llabel(f"{ele_list[1]}")
	ax.set_rlabel(f"{ele_list[2]}")
	if equi:
		plt.savefig(f"{folder_path}/{ele_list}/{temp}_equi.png")
	else:
		plt.savefig(f"{folder_path}/{ele_list}/{temp}.png")
