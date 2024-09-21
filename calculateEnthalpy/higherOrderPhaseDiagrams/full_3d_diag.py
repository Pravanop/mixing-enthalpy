import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpltern
from tqdm import tqdm
from calculateEnthalpy.helper_functions.grid_code import create_mol_grid
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram



correction = True
equi = False
im_flag = True
isShow = False
cmap = plt.get_cmap('viridis')

if correction:
	binary_file_path = "../../calcEnthalpy_old/new_phase_diagram/bokas_omegas_processed.json"
else:
	binary_file_path = "../data/output_data/bokasCorrected_bcc_1/all_lattices_binaries.json"

end_member_path = "../../calcEnthalpy_old/new_phase_diagram/bokas_end_members_dict.json"

pD = phaseDiagram(
	processed_binary_file_path=binary_file_path,
	end_member_file_path=end_member_path,
	grid_size=15,
	im_flag=im_flag,
	correction=correction,
	equi_flag=equi)

composition = ['Cr', 'V', 'Ti']
composition.sort()
temp_grid = np.linspace(0, 3500, 20).astype(int)
phase_diag_dict = pD.make_PD_comp_temp(composition=composition, temp_grid=temp_grid)

temps = list(phase_diag_dict.keys())
grid_size = 20
N_i = 100
skip = 1
sigma = 1
levels = 5

mol_grid = create_mol_grid(3, grid_size)

t, l, r = mol_grid[:, 0], mol_grid[:, 1], mol_grid[:, 2]

temp_stable = {}
for temp in tqdm(temps, desc="Creating phase diagrams"):
	stables = []
	for idx, mol in enumerate(mol_grid):
		enthalpy, entropy, mol_ratio = pD.find_enthalpy_entropy_composition(composition=composition,
																   mol_ratio=mol,
																	lattice = 'min',
																	temperature=temp)
		is_stable = pD.check_stability(mol_ratio=mol_ratio,
									   temp=temp,
									   conv_hull=phase_diag_dict[temp],
									   entropy=entropy,
									   mix_enthalpy=enthalpy,)
		
		if is_stable is not None:
			if 'Ti' not in mol_ratio.keys():
				print(mol_ratio, is_stable[1], temp)
			if np.isclose(np.round(is_stable[1],4), 0.0, atol=1e-3):
				stable = 0
			else:
				stable = 1
		else:
			if 'Ti' not in mol_ratio.keys():
				print(mol_ratio, is_stable, temp)
			stable = 0
		
		if 'Ti' not in mol_ratio.keys():
			print(mol_ratio, stable, temp)
		stables.append(stable)
		
	# break
	
	temp_stable[temp] = stables

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i, stable in temp_stable.items():
	ax.scatter(mol_grid[:, 0], mol_grid[:, 1], c = stable, cmap = 'coolwarm', zs=i, zdir='z', label=f'Set {i + 1}')

ax.grid(False)
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.view_init(elev=30, azim=270)
plt.show()
