import pandas as pd
import numpy as np

from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram

file_path = "/Users/pravanomprakash/Library/CloudStorage/Box-Box/HEA_alloys/sheets/HEA candidates enthalpy data.xlsx"
df_master = pd.read_excel(file_path, sheet_name=0)
index = np.array(df_master.iloc[:, -1])
names = np.array(df_master.iloc[:, 0])
total = np.vstack((index, names)).T
master_idx = []
for i in total:
	if not np.isnan(i[0]):
		master_idx.append(i)

master_idx = np.array(master_idx)
processed_file_path = "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/data/output_data/Pravan_bcc_4/all_lattices_binaries.json"
pD = phaseDiagram(
	processed_file_path=processed_file_path,
grid_size=20)

# composition = ['Cr', 'Fe', 'Mn']
equi = True


# #Find Miscibility temperature
calphad = []
model = []
for ind, k in enumerate(master_idx):
	df_ind = pd.read_excel(file_path, sheet_name=int(k[0]))
	phase_info = df_ind.loc[:, 'H':]
	temps = df_ind['T'].to_numpy()
	df_fraction = df_ind.filter(regex=("f.*"))
	cols = df_fraction.columns.to_numpy()
	for idx2, j in enumerate(cols):
		cols[idx2] = j[3:-1]

	results = []
	misc_temp = []
	for idx, i in enumerate(temps):
		values = df_fraction.iloc[idx].to_numpy()
		tieline = np.setxor1d(np.nonzero(values), np.argwhere(np.isnan(values)).flatten())
		phases = cols[tieline]
		if phases.size == 1:
			if 'Liquid' not in phases:
				misc_temp.append(i)

	composition = k[1].split('-')
	if equi:
		mol_ratio = [1 / len(composition)] * len(composition)
		# results.append(f'{i}K : {"-".join(phases)}')
	misc_T = pD.find_misc_temperature(composition=composition, mol_ratio=mol_ratio)
	misc_temp = max(misc_temp) if misc_temp else None
	calphad.append(misc_temp)
	model.append(misc_T)

df_results = pd.DataFrame([master_idx[:, 1], calphad, model]).T
df_results.columns = ['Alloy', 'Calphad', 'Model']
print(df_results)
df_results.to_csv("Pravan_bcc_4_calphad.csv", index=False)