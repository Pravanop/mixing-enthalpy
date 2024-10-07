import pandas as pd
from far_heaa.io.json_handler import JSONHandler
import numpy as np
import os
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt




class SubRegularExtraction:
	"""
	TODO Create a more stringent database for sub-regular.
	TODO Calculate the end member energies and change the end_member dict accordingly.
	TODO Create the input csv with sorted binary data. Order Matters!
	"""

	def __init__(self):
		self.df = pd.read_csv('../database/Outputs_BCC_runs.csv')
		self.lattice = 'BCC'
		self.end_member_data = JSONHandler.load_json(
			folder_path='/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calcEnthalpy_old/new_phase_diagram/',
			file_name='bokas_omegas')
		# self.temporary_corrections()
	@staticmethod
	def sub_regular_model(x_i, omega1, omega2):
		return x_i * (1 - x_i) * (omega1 * x_i + omega2 * (1 - x_i))

	def temporary_corrections(self):
		#for BCC unit end_member data is made by us
		self.end_member_data = self.end_member_data['elements'][self.lattice]
		self.end_member_data['Cr'] = -9.52
		self.end_member_data['W'] = -12.96
		self.end_member_data['Ta'] = -11.82
		self.end_member_data['V'] = -8.97
# alloys = sub_regular_binary_bcc['Unnamed: 0']

	def store_partial_sub_regular(self):
		h_mix = []
		mol_ratios = self.df.columns[1:].to_numpy().astype(float)
		mol_ratios = np.append(0, mol_ratios)
		mol_ratios = np.append(mol_ratios, 1)
		popts = {}
		x = np.linspace(0, 1, 100)
		for idx, row in self.df.iterrows():
			h_dft = row[1:].to_numpy().astype(float)
			# h_dft = h_dft / 24

			end_members = row[0].split('-')
			end_members.reverse()
			end_member_energies = np.array([self.end_member_data[end_member] for end_member in end_members])
			# print()
			h_dft = np.pad(h_dft, (1, 1), 'constant')
			is_nan = np.isnan(h_dft)
			y_data = h_dft.copy()
			y_data = y_data[~is_nan]
			x_data = mol_ratios.copy()
			x_data = x_data[~is_nan]
			sorted_endmembers = sorted(end_members)
			y_data = y_data - x_data * end_member_energies[0] - (1 - x_data) * end_member_energies[1]
			y_data[0] = 0
			y_data[-1] = 0
			h_mix.append(y_data)
			popt, pcov = curve_fit(self.sub_regular_model, x_data, y_data)
			y = self.sub_regular_model(x_data, *popt)
			assert np.log10(np.linalg.cond(pcov)) < 10
			popt = np.round(popt, 8)
			popts['-'.join(sorted_endmembers)] = popt

		binary_dict = JSONHandler.load_json(folder_path='../database/', file_name='bokas_omegas_processed')

		for key, value in binary_dict.items():
			if key in popts:
				print(key)
				binary_dict[key][self.lattice] = list(popts[key])

		JSONHandler.save_json(data=binary_dict, folder_path='../database/', file_name='bokas_omegas_processed_subregular')
