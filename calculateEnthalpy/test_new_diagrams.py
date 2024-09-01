import itertools

import pandas as pd

from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram



def avg_T_melt(composition):
	meltT = pd.read_csv("../data/input_data/PubChemElements_all.csv").to_numpy()
	meltT = dict(zip(meltT[:, 0], meltT[:, 1]))
	tm = [meltT[i] for i in composition]
	return int(sum(tm) / len(tm))

# Cr,Ta,Ti,V,W
compositions = list(itertools.combinations(['Cr', 'Fe', 'Mn', 'Hf', 'V', 'W','Ti', 'Ta'], r = 4))
compositions = ['-'.join(sorted(i)) for i in compositions]
with_bias = []
without_bias = []
avg_T = []
for composition in compositions:
	print(composition)
	composition = composition.split('-')
	print("Average T_m:", avg_T_melt(composition))

	equi = True
	if equi:
		mol_ratio = [1 / len(composition)] * len(composition)
	avg_T.append(avg_T_melt(composition))
	for correction in [True, False]:

		if correction:
			binary_file_path = "new_phase_diagram/bokas_omegas_processed.json"
		else:
			binary_file_path = "../data/output_data/bokasCorrected_bcc_1/all_lattices_binaries.json"

		end_member_path = "new_phase_diagram/bokas_end_members_dict.json"

		pD = phaseDiagram(
			processed_binary_file_path=binary_file_path,
			end_member_file_path=end_member_path,
			grid_size=10,
			im_flag=True,
			correction=correction,
			equi_flag=False)

		# print("Biased:", correction)
		for lattice in ['BCC']:
			misc_T = pD.find_misc_temperature(composition=composition, lattice=lattice, mol_ratio=mol_ratio)
			if correction:
				with_bias.append(misc_T)
			else:
				without_bias.append(misc_T)
			print(f"T_misc: (bias: {correction})", misc_T)

compositions = pd.Series(compositions, name='compositions').T
with_bias = pd.Series(with_bias, name='with_bias').T
without_bias = pd.Series(without_bias, name='without_bias').T
avg_T = pd.Series(avg_T, name='avg_T').T
print(with_bias)
df = pd.concat([compositions, with_bias, without_bias, avg_T], axis=1)
# df.columns = ['composition', 'with_bias','avg_T', 'without_bias']

print(df.head())
df.to_csv("../data/output_data/quarternary_predictions_bokas.csv")