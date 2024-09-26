import pandas as pd
from far_heaa.io.json_handler import JSONHandler
import numpy as np
import os
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

df = pd.read_csv('../database/binary_bcc.csv')
lattice = 'BCC'
end_member_data = JSONHandler.load_json(
	folder_path='/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calcEnthalpy_old/new_phase_diagram/',
	file_name='bokas_omegas')


def sub_regular_model(x_i, omega1, omega2):
	return x_i * (1 - x_i) * (omega1 * x_i + omega2 * (1 - x_i))


end_member_data = end_member_data['elements'][lattice]
end_member_data['Cr'] = -9.52
end_member_data['W'] = -12.96
end_member_data['Ta'] = -11.82
end_member_data['V'] = -8.97
# alloys = sub_regular_binary_bcc['Unnamed: 0']

h_mix = []
mol_ratios = df.columns[1:].to_numpy().astype(float)
mol_ratios = np.append(0, mol_ratios)
mol_ratios = np.append(mol_ratios, 1)
inverse_mol_ratios = 1 - mol_ratios
popts = {}
x = np.linspace(0, 1, 100)
for idx, row in df.iterrows():
	h_dft = row[1:].to_numpy().astype(float)
	h_dft = h_dft / 24
	
	end_members = row[0].split('-')
	end_members.reverse()
	end_member_energies = np.array([end_member_data[end_member] for end_member in end_members])
	print(end_member_energies, end_members)
	# print()
	h_dft = np.pad(h_dft, (1, 1), 'constant')
	is_nan = np.isnan(h_dft)
	y_data = h_dft.copy()
	y_data = y_data[~is_nan]
	x_data = mol_ratios.copy()
	x_data = x_data[~is_nan]
	sorted_endmembers = sorted(end_members)
	print(y_data)
	y_data = y_data - x_data * end_member_energies[0] - (1 - x_data) * end_member_energies[1]
	y_data[0] = 0
	y_data[-1] = 0
	h_mix.append(y_data)
	# print(end_member_energies)
	popt, pcov = curve_fit(sub_regular_model, x_data, y_data)
	y = sub_regular_model(x_data, *popt)
	print(np.linalg.norm(y_data - y))
	assert np.log10(np.linalg.cond(pcov)) < 10
	popt = np.round(popt, 8)
	popts['-'.join(sorted_endmembers)] = popt

# y = sub_regular_model(x, *popts['Cr-W'])
# # plt.plot(y, 'o', label='data')
# # plt.show()
binary_dict = JSONHandler.load_json(folder_path='../database/', file_name='bokas_omegas_processed')

for key, value in binary_dict.items():
	if key in popts:
		print(key)
		binary_dict[key][lattice] = list(popts[key])

JSONHandler.save_json(data=binary_dict, folder_path='../database/', file_name='bokas_omegas_processed_subregular')
