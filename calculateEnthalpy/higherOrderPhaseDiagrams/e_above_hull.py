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
# composition.sort()
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
e_above_hull = []
mol = [0.6667, 0.3333, 0.0000]


for temp in tqdm(temps, desc="Creating phase diagrams"):
	enthalpy, entropy, mol_ratio = pD.find_enthalpy_entropy_composition(composition=composition,
																		mol_ratio=mol,
																		lattice='min',
																		temperature=temp)
	print(mol_ratio, enthalpy, entropy, enthalpy - temp*entropy)
	is_stable = pD.check_stability(mol_ratio=mol_ratio,
								   temp=temp,
								   conv_hull=phase_diag_dict[temp],
								   entropy=entropy,
								   mix_enthalpy=enthalpy,)
		
	if is_stable is not None:
		e_above_hull.append(np.round(is_stable[1] ,4))
	else:
		e_above_hull.append(0)
		
plt.scatter(temps, e_above_hull)
plt.ylabel('E above hull (eV/atom)')
plt.xlabel('Temperature (K)')
plt.show()
