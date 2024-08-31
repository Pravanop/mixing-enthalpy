import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
matplotlib.rcParams.update({'font.size': 14})
omegas_dict = json.load(open('./bokas_omegas_processed.json'))
end_members_dict = json.load(open('./bokas_end_members_dict.json'))

binary = 'Cr-Ti'
lattices = ['BCC', 'HCP']
omega_info = omegas_dict[binary]
end_members = sorted(binary.split('-'))
end_members_info = [end_members_dict[i] for i in end_members]

print(omega_info)
# print(end_members_info, omega_info)

def calculate_binary_mixing_enthalpy(binary_list, omega_info, end_members_info):
	mix_enthalpy = {}
	mix_enthalpy_original = {}
	x = np.linspace(0, 1, 100)

	for lattice in lattices:
		reg_solution = 0
		reg_solution = omega_info[lattice] * x * (1 - x)  #regular solution
		# print(reg_solution[1])
		mix_enthalpy_original[lattice] = reg_solution
		#biasing
		reg_solution2 = reg_solution +  end_members_info[0][lattice]
		#adding slope
		reg_solution2 += (end_members_info[1][lattice] - end_members_info[0][lattice])*x
		mix_enthalpy[lattice] = reg_solution2

	return mix_enthalpy

lattice_colors = dict(zip(lattices, ['r','b']))
def hmix_plotter(mixing_enthalpy_dict):
	for key, enthalpy in mixing_enthalpy_dict.items():
		plt.plot(x, enthalpy, label = key, color = lattice_colors[key])

		plt.xticks([0, 1], labels=end_members)
		plt.legend()
		plt.xlim([0,1])
		plt.ylim([-0.1, 0.1])
		plt.xlabel('$x$')
		plt.ylabel('$H_{mix}$ (eV/atom)')
		plt.savefig('mix_enthalpy.png', dpi = 300)

# hmix_plotter(mix_enthalpy_original)
kb = 8.61773e-5
config_entropy = -kb*(x*np.log(x) + (1-x)*np.log(1-x))

T_grid = np.linspace(0, 3200, 10).astype(int).reshape(1,-1).T
def gibbs(enthlapy, temperature, entropy):

	# print(temperature*entropy)
	return enthlapy - temperature*entropy

gibbs_dict = {}
for key, value in mix_enthalpy_original.items():
	gibbs_lattice = gibbs(value, T_grid, config_entropy)
	gibbs_dict[key] = gibbs_lattice

def gibbs_plotter(gibbs_dict):
	for key, value in gibbs_dict.items():
		plt.axhline(y=0, color='#BBBBBB', linestyle='-')
		for idx, i in enumerate(value):
			plt.plot(x, i , label = T_grid[idx][0], color = lattice_colors[key])

	plt.xlim([0,1])
	# plt.ylim([-1, 1])
	plt.xticks([0,1], labels = end_members)
	plt.xlabel('$x$')
	plt.ylabel('$G$ (eV/atom)')
	plt.subplots_adjust(right=0.9, left=0.2)
	patches = []
	for key, value in lattice_colors.items():
		patches.append(mpatches.Patch(color=value, label=key))
	plt.legend(handles=patches)
	plt.savefig('gibbs.png', dpi = 300)

gibbs_plotter(gibbs_dict)