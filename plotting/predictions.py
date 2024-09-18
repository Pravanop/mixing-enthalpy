import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm

from calculateEnthalpy.helper_functions.grid_code import create_multinary, create_mol_grid
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram
import matplotlib

matplotlib.rcParams['font.size'] = 24
correction = True
equi = False

if correction:
	binary_file_path = "../calcEnthalpy_old/new_phase_diagram/bokas_omegas_processed.json"
else:
	binary_file_path = "../../data/output_data/bokasCorrected_bcc_1/all_lattices_binaries.json"

end_member_path = "../calcEnthalpy_old/new_phase_diagram/bokas_end_members_dict.json"
pD = phaseDiagram(
	processed_binary_file_path=binary_file_path,
	end_member_file_path=end_member_path,
	grid_size=10,
	im_flag=True,
	correction=correction,
	equi_flag=equi)


dim = 5

compositions = create_multinary(element_list=['Cr', 'Ta', 'Ti', 'V', 'W', 'Hf', 'Fe', 'Mn'], no_comb=[dim])
compositions = list(compositions.values())[0]
acg_tm = []
temp_list = []
alloys = []
for composition in tqdm(compositions, desc='Calculating Misc_T'):
	composition = composition.split('-')
	n_alloy = len(composition)
	all_combs = create_multinary(element_list=composition, no_comb=list(range(2, n_alloy + 1)))
	im_list = []
	for dimensionality, alloy_list in all_combs.items():
		if pD.im_flag:
			im_list += pD.get_intermetallic(alloy_list)

	mol = [1/len(composition)]*len(composition)
	lattice = "HCP"
	misc_T = pD.find_misc_temperature(composition=composition, lattice=lattice, mol_ratio=mol, batch_tag=True, im=im_list, phase_flag=False)
	if isinstance(misc_T, float):
		temp_list.append(misc_T)
	else:
		temp_list.append(misc_T[0])
	alloys.append('-'.join(composition))
	acg_tm.append(pD.temp_grid[-1])

df = pd.DataFrame([alloys, temp_list, acg_tm])
df = df.T
df.columns = ["Alloys", 'Temp', 'Avg_Tm']
df_proc = df.loc[df['Temp'] < df['Avg_Tm']]
df.to_csv(f"./misc_T_{dim}_equi_{lattice}.csv", index=False)
print(df_proc.shape)
df_proc = df_proc.sort_values(['Avg_Tm'])
fig, ax = plt.subplots(figsize = (18, 8))
ax.bar(x = df_proc['Alloys'], height = df_proc['Avg_Tm'], color = '#4393C3', edgecolor = 'black', width = 0.5, align='center', zorder=2)
ax.bar(x = df_proc['Alloys'], height = df_proc['Temp'], color = '#D6604D', edgecolor = 'black', width = 0.5, align='center', zorder=2)
ax.axhline(y = 2000, color = 'black', linestyle = '--', linewidth = 3.5, alpha = 0.8, zorder=0, label = '_nolegend_')
ax.legend(['Average T$_{m}$', 'Miscible Temperature', ])
# ax.set_xlabel('Alloys')
ax.set_ylabel('T (K)')
ax.tick_params(axis='x', labelrotation=90)
plt.subplots_adjust(bottom=0.3, top = 0.99, left=0.1, right=0.99)
if equi:
	plt.savefig(f'../plots/misc_T_{dim}_equi_{lattice}.png')
else:
	plt.savefig(f'../plots/misc_T_{dim}_{lattice}.png')
