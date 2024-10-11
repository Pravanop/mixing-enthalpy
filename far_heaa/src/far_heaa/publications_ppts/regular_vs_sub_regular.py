import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from far_heaa.io.json_handler import JSONHandler
import matplotlib.pyplot as plt
import seaborn as sns

end_member_data = {
	
	"Cr": "-9.52005384",
	"Fe": "-8.25682767",
	"Hf": "-9.74533769",
	"Mo": "-10.93282547",
	"Nb": "-10.10322176",
	"Ta": "-11.81975660",
	"Ti": "-7.69671714",
	"V": "-8.96529006",
	"W": "-12.95591411",
	"Zr": "-8.43952352"
	
}


def sub_regular_model(x_i, omega1, omega2):
	return x_i * (1 - x_i) * (omega1 * x_i + omega2 * (1 - x_i))

def regular_model(x_i, omega):
	return x_i * (1 - x_i) * omega

df = pd.read_csv('../database/Outputs_BCC_runs.csv')
mol_ratios = df.columns[1:].to_numpy().astype(float)
mol_ratios = np.append(0, mol_ratios)
mol_ratios = np.append(mol_ratios, 1)
popts = {}
h_mix = []
binary_dict = JSONHandler.load_json(folder_path='../database/', file_name='bokas_omegas_processed')
regular_popts = {}
for idx, row in df.iterrows():
	h_dft = row[1:].to_numpy().astype(float)
	# reg_mix = h_mix[2]
	end_members = row[0].split('-')
	end_members.reverse()
	end_member_energies = np.array([float(end_member_data[end_member]) for end_member in end_members])
	h_dft = np.pad(h_dft, (1, 1), 'constant')
	is_nan = np.isnan(h_dft)
	y_data = h_dft.copy()
	y_data = y_data[~is_nan]
	x_data = mol_ratios.copy()
	x_data = x_data[~is_nan]/100
	y_data = y_data - x_data * end_member_energies[0] - (1 - x_data) * end_member_energies[1]
	y_data[0] = 0
	y_data[-1] = 0
	h_mix.append(y_data)
	
	popt, pcov = curve_fit(sub_regular_model, x_data, y_data)
	assert np.log10(np.linalg.cond(pcov)) < 10
	popt = np.round(popt, 8)
	popts[row[0]] = popt

	try:
		reg_mix = y_data[np.where(x_data == 0.5)[0][0]]
	except:
		reg_mix = sub_regular_model(0.5, *popt)
		print(0)

	regular_popts[row[0]] = reg_mix*4


x = np.linspace(0, 1, 100)
x_real = [0, 0.125, 0.25, 0.5, 0.75, 0.875, 1]
asymmetric_list = []
diff_omega = []
for idx, key in enumerate(popts.keys()):
	sub_reg_y = sub_regular_model(x, *popts[key])
	# plt.plot(x, sub_reg_y, label=key)
	# plt.scatter(x_real, h_mix[idx])
	# plt.title(key)
	# plt.show()
	reg_y = regular_popts[key]*x*(1-x)
	diff_omega.append(popts[key][1] - popts[key][0])
	extremum_reg = np.max(np.abs(reg_y))
	extremum_reg_x = x[np.argmax(np.abs(reg_y))]
	extremum_sub_reg = np.max(np.abs(sub_reg_y))
	extremum_sub_reg_x = x[np.argmax(np.abs(sub_reg_y))]
	asymmetric_dist = np.sqrt((extremum_reg_x - extremum_sub_reg_x)**2 + (extremum_reg - extremum_sub_reg)**2)
	asymmetric_list.append(np.round(asymmetric_dist,2))

asymmetric_list = np.array(asymmetric_list)
diff_omega = np.array(diff_omega)
#normalize asymmetric distances
asymmetric_list = np.round((asymmetric_list - np.min(asymmetric_list))/(np.max(asymmetric_list) - np.min(asymmetric_list)),2)
# diff_omega = np.round((diff_omega - np.min(diff_omega))/(np.max(diff_omega) - np.min(diff_omega)),2)

element_list = np.unique(np.array([key.split('-') for key in popts.keys()]).flatten())
print(element_list)
element_list = np.array(['Mo','Nb','Ta','V','W','Zr','Ti'])
#make a heatmap of the asymmetric distances
asymmetric_matrix = np.zeros((len(element_list), len(element_list)))
for idx, key in enumerate(popts.keys()):
	elements = key.split('-')
	element_idx = [np.where(element_list == element)[0][0] for element in elements]
	asymmetric_matrix[element_idx[0], element_idx[1]] = diff_omega[idx]
	# asymmetric_matrix[element_idx[0], element_idx[1]] = asymmetric_list[idx]
	asymmetric_matrix[element_idx[1], element_idx[0]] = diff_omega[idx]
	# asymmetric_matrix[element_idx[1], element_idx[0]] = asymmetric_list[idx]

sns.heatmap(asymmetric_matrix, annot=True,
			cmap='plasma',
			xticklabels=element_list,
			yticklabels=element_list,
			cbar_kws={'label': 'Asymmetric Distance'},)


plt.savefig('../plots/asymmetric_distances_omega.png', dpi=100)

	