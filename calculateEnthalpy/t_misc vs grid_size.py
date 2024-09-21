from matplotlib import pyplot as plt

from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
correction = False
if correction:
	binary_file_path = "../calcEnthalpy_old/new_phase_diagram/bokas_omegas_processed.json"
else:
	binary_file_path = "data/output_data/bokasCorrected_bcc_1/all_lattices_binaries.json"

end_member_path = "../calcEnthalpy_old/new_phase_diagram/bokas_end_members_dict.json"

misc_T = []
composition = ['Cr', 'W', 'Fe', 'Ti']
mol_ratio = [1/len(composition)]*len(composition)

for grid_size in range(3,31,2):
	print(grid_size)
	pD = phaseDiagram(
		processed_binary_file_path=binary_file_path,
		end_member_file_path=end_member_path,
		grid_size=grid_size,
		im_flag=False,
		correction=correction,
		equi_flag=False)
	misc_T.append(pD.find_misc_temperature(composition=composition,
										   mol_ratio=mol_ratio,
										   lattice='BCC'))

plt.plot(range(3,31,2), misc_T, marker='o', color='black')
plt.xlim(3,30)
plt.ylabel('Temperature (K)')
plt.xlabel('Grid Density')
plt.title('-'.join(composition))
plt.subplots_adjust(left=0.15)
plt.savefig('misc_temperature_vs_grid_size.png', dpi = 300)