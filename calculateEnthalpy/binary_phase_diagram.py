import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm

from calculateEnthalpy.helper_functions.grid_code import create_mol_grid, create_multinary
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram
import matplotlib
from scipy.interpolate import CubicSpline, splrep, BSpline, UnivariateSpline

matplotlib.rcParams.update({'font.size': 18})
correction = True
if correction:
	binary_file_path = "new_phase_diagram/bokas_omegas_processed.json"
else:
	binary_file_path = "../data/output_data/bokasCorrected_bcc_1/all_lattices_binaries.json"

end_member_path = "new_phase_diagram/bokas_end_members_dict.json"

misc_T = []
composition = ['Cr', 'Ta']
composition = sorted(composition)
mol_grid = create_mol_grid(n=2, grid_size=30)
pD = phaseDiagram(
	processed_binary_file_path=binary_file_path,
	end_member_file_path=end_member_path,
	grid_size=30,
	im_flag=True,
	correction=correction,
	equi_flag=False)

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

for i in tqdm(mol_grid):

	misc_T = pD.find_misc_temperature(composition=composition, mol_ratio=i, lattice='BCC', phase_flag=True, batch_tag=True, im=im_list)
	if isinstance(misc_T, float):
		misc_T_list.append(misc_T)
	else:
		misc_T_list.append(5000)

cs = UnivariateSpline(mol_grid[:, 0], misc_T_list)
cs.set_smoothing_factor(2)
x_range = np.arange(0, 1.01, 0.01)
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_facecolor('#004488')
ax.plot(x_range, cs(x_range), label='Cubic Spline', color='black', linewidth=4)
ax.fill(x_range, cs(x_range), color='#DDAA33')

# plt.scatter(x=mol_grid[:, 0], y=misc_T_list)
ax.set_xlim(0,1)
ax.set_xticks([0,1], sorted(composition, reverse=True))
ax.set_ylabel('T (K)')
ax.set_xlabel('x')
ax.set_yticks(np.linspace(0,3000, 6))
ax.set_ylim(0,3000)
ax.set_title('-'.join(sorted(composition, reverse=True)) + " Phase Diagram")
plt.subplots_adjust(left = 0.15, right=0.9)
plt.savefig(f'binary_diagram_{"-".join(sorted(composition, reverse=True))}.png', dpi=300)
