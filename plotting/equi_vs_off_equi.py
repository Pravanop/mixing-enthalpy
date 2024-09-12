import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm

from calculateEnthalpy.helper_functions.grid_code import create_multinary, create_mol_grid
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram
import matplotlib

matplotlib.rcParams['font.size'] = 16
correction = True
# equi = False

if correction:
	binary_file_path = "../calcEnthalpy_old/new_phase_diagram/bokas_omegas_processed.json"
else:
	binary_file_path = "../../data/output_data/bokasCorrected_bcc_1/all_lattices_binaries.json"

end_member_path = "../calcEnthalpy_old/new_phase_diagram/bokas_end_members_dict.json"

equimolar_temps = []
non_equimolar_temps = []
dims = []
for equi in [False, True]:
	pD = phaseDiagram(
		processed_binary_file_path=binary_file_path,
		end_member_file_path=end_member_path,
		grid_size=10,
		im_flag=True,
		correction=correction,
		equi_flag=equi)

	for dim in [2, 3]:

		compositions = create_multinary(element_list=['Cr', 'Ta', 'Ti', 'V', 'W', 'Hf', 'Fe', 'Mn'], no_comb=[dim])
		compositions = list(compositions.values())[0]

		for composition in tqdm(compositions, desc='Calculating Misc_T'):
			if equi:
				dims.append(dim)
			composition = composition.split('-')
			n_alloy = len(composition)
			all_combs = create_multinary(element_list=composition, no_comb=list(range(2, n_alloy + 1)))
			im_list = []

			for dimensionality, alloy_list in all_combs.items():
				if pD.im_flag:
					im_list += pD.get_intermetallic(alloy_list)

			mol = [1 / len(composition)] * len(composition)
			misc_T = pD.find_misc_temperature(composition=composition, lattice="min", mol_ratio=mol, batch_tag=True,
											  im=im_list, phase_flag=False)
			if isinstance(misc_T, float):
				if equi:
					equimolar_temps.append(misc_T)
				else:
					non_equimolar_temps.append(misc_T)
			else:
				if equi:
					equimolar_temps.append(misc_T[0])
				else:
					non_equimolar_temps.append(misc_T[0])

fig, ax = plt.subplots()
print(len(equimolar_temps), len(non_equimolar_temps), len(dims))
ax.scatter(equimolar_temps, non_equimolar_temps, c=dims, zorder=1)
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.plot([0, 3500], [0, 3500], linestyle='--', c='black', linewidth=1, label='_nolegend_', zorder=0)
ax.set_xlim([0, 3500])
ax.set_ylim([0, 3500])
ax.set_xlabel('T$_{misc}$ [equi] (K)')
ax.set_ylabel('T$_{misc}$ [off_equi] (K)')
plt.subplots_adjust(left=0.2, bottom=0.15, right=0.95)
plt.savefig('../plots/equi_vs_off_equi.png', dpi=300)
