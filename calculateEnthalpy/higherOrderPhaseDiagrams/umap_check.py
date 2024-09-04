import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm
import mpltern
from calculateEnthalpy.helper_functions.grid_code import create_multinary, create_mol_grid
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram

composition = ['Cr', 'V', 'Ti']

tm = []
n_alloy = len(composition)
all_combs = create_multinary(element_list=composition, no_comb=list(range(2, n_alloy + 1)))
im_list = []

mol_grid = create_mol_grid(3, 20)
mol_grid = np.array(mol_grid)

correction = True
equi = False

if correction:
	binary_file_path = "../new_phase_diagram/bokas_omegas_processed.json"
else:
	binary_file_path = "../../data/output_data/bokasCorrected_bcc_1/all_lattices_binaries.json"

end_member_path = "../new_phase_diagram/bokas_end_members_dict.json"

pD = phaseDiagram(
	processed_binary_file_path=binary_file_path,
	end_member_file_path=end_member_path,
	grid_size=20,
	im_flag=True,
	correction=correction,
	equi_flag=equi)

for dimensionality, alloy_list in all_combs.items():
	if pD.im_flag:
		im_list += pD.get_intermetallic(alloy_list)
for ind, mol in enumerate(tqdm(mol_grid, desc = "Iterating through phase space")):
	# mol_ratio = [i[str(j)] for j in range(1, n+1)]
	for lattice in ['min']:
		misc_T = pD.find_misc_temperature(composition=composition, lattice=lattice, mol_ratio=mol, batch_tag=True, im=im_list)
		if isinstance(misc_T, float):
			tm.append(misc_T)
		else:
			tm.append(-1)


t, l, r = mol_grid[:, 0], mol_grid[:, 1], mol_grid[:, 2]
fig = plt.figure()

ax = fig.add_subplot(projection="ternary")
ax.grid()
cax = ax.inset_axes([1.03, 0.1, 0.05, 0.9], transform=ax.transAxes)
# pc = ax.scatter(t, l, r, c = tm, cmap = "coolwarm")
# data = np.concat([t, l, r, ])
ps = ax.tricontour(t,l,r,tm, levels=10, cmap = "coolwarm")
colorbar = fig.colorbar(ps, cax=cax)
colorbar.set_label('T_misc', rotation=270, va='baseline')
ax.set_tlabel(f"{composition[0]}")
ax.set_llabel(f"{composition[1]}")
ax.set_rlabel(f"{composition[2]}")

plt.savefig(f"../../plots/phase_diagrams/{'-'.join(composition)}")
# plt.show()
