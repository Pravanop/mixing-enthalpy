import pickle

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import accuracy_score, confusion_matrix, r2_score
from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.io.json_handler import JSONHandler
from far_heaa.high_throughput.realignment_maths import get_realignment_energy
import seaborn as sns

from far_heaa.math_operations.thermo_calculations import ThermoMaths

system = 3
element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf']
file = JSONHandler.load_json(folder_path='/Users/pravanomprakash/Library/CloudStorage/Box-Box/HEA_alloys/done_data/', file_name=f'{system}_add_ele_paths_total_10_BCC_wo_correction')

data = JSONHandler.load_json(folder_path='../database', file_name='bokas_omegas_processed')
mol_grid_size = 5
x = np.linspace(0, 1, mol_grid_size)

results_dict = {}
miscible_keys = []
never_miscible_keys = []
possible_miscible_keys = []
tm = ThermoMaths()
throw_away = []
error_count = 0

h_alloy_list = []
h_comp_list = []
h_r_list = []
temperature_list = []
im_formation_list = []
reduction = []
with open('../database/intermetallic_database.pickle', 'rb') as handle:
	im_list = pickle.load(handle)

for key, value in file.items():
	for add_ele, temp_list in value.items():
		composition = key.split('-')
		if add_ele in composition:
			composition.remove(add_ele)
		temp_array = np.array(temp_list)
		temp_array[temp_array == -1000.0] = np.nan
		
		if np.isnan(temp_array[0]) or np.isnan(temp_array[-1]):
			continue
		
		h_alloy = 0
		im_formation = []
		
		for j in composition:
			h_alloy += data['-'.join(sorted([j, add_ele]))]['BCC']
			im_formation += [k.energy / k.composition.num_atoms for k in im_list['-'.join(sorted([j, add_ele]))] if
							 '-'.join(sorted([j, add_ele])) in im_list]
		
		h_alloy_list.append(h_alloy)
		h_r = get_realignment_energy(element_list=composition + [add_ele], lattice='BCC', data=data)
		h_r_list.append(h_r)
		h_comp = 0
		for j in list(MultinaryCombinations.create_multinary(composition, no_comb=[2]).values())[0]:
			h_comp += data[j]['BCC']
			# im_formation += [k.energy / k.composition.num_atoms for k in im_list[j] if j in im_list]
		if '-'.join(sorted(composition + [add_ele])) in im_list:
			im_formation += [k.energy / k.composition.num_atoms for k in
							 im_list['-'.join(sorted(composition + [add_ele]))]]
        
		if not im_formation:
			im_formation = 0
		h_comp_list.append(h_comp)
		im_formation_list.append(min(im_formation))
		
		alloy_temp = tm.avg_T_melt(composition=composition, mol_ratio=[1 / system] * system)
		element_temp = tm.avg_T_melt(composition=composition + [add_ele], mol_ratio=[1 / (system + 1)] * (system + 1))
		temperature_list.append(float(alloy_temp - element_temp))
		if np.isnan(temp_array[-1]):
			reduction.append(-2000 - temp_array[0])
		elif np.isnan(temp_array[0]):
			reduction.append(temp_array[-1] + 3500)
		else:
			reduction.append(temp_array[-1] - temp_array[0])

reduction_flag = [0 if i < 0 else 1 for i in reduction]
# reduction_flag = np.array(reduction)
# reduction_flag = (reduction_flag - min(reduction_flag))/(max(reduction_flag) - min(reduction_flag))
df = pd.DataFrame([h_alloy_list, h_comp_list,h_r_list, temperature_list, im_formation_list, reduction_flag]).T
# df = pd.DataFrame([h_alloy_list, h_comp_list, temperature_list, im_formation_list, reduction_flag]).T
total_columns = ['Ha', 'Hc', 'Hr', 'delta meltT', 'im', 'Miscibility']
# total_columns = ['Ha', 'Hc', 'h_r', 'im', 'Miscibility']
df.columns = total_columns

fig, axs = plt.subplots(len(total_columns[:-1]), 1, sharey=True)

# Plotting each feature in a separate subplot
features = total_columns[:-1]
for i, feature in enumerate(features):
	sns.histplot(df[feature], bins=15, color='#BB5566', edgecolor='black', stat='percent', ax=axs[i])
	axs[i].set_xlabel(feature)
	axs[i].set_ylabel('Frequency')

plt.subplots_adjust(wspace=0, hspace=0)
# Adjust layout for better spacing
plt.tight_layout()
plt.show()
# plt.savefig(f'../plots/high_throughput/descriptor_featdist_{system}.png', dpi=150)

pairplot = sns.pairplot(df, hue='Miscibility')
legend = pairplot._legend
legend.set_title("Miscibility")  # Set a custom title
legend.get_texts()[0].set_text("Increases")  # Change the text for the first group
legend.get_texts()[1].set_text("Decreases")
# plt.show()
# plt.savefig(f'../plots/high_throughput/descriptor_pairplot_{system}.png', dpi=150)
#
fig, ax = plt.subplots()
sns.heatmap(df.corr(), annot=True,
			xticklabels=['Ha', 'Hc', 'delta meltT', 'im' 'Miscibility'],
			yticklabels=['Ha', 'Hc', 'delta meltT', 'im' 'Miscibility'],
			ax=ax)
print(df.info())
plt.show()
# plt.savefig(f'../plots/high_throughput/descriptor_corrmat_{system}.png', dpi=150)

model = LogisticRegression()
# model = LinearRegression()
X = df[features].to_numpy()
y = df['Miscibility'].to_numpy()
print(X.shape, y.shape)
model.fit(X, y)
coef = model.coef_[0]
intercept = model.intercept_[0]

# Example prediction
y_pred = model.predict(X)
accuracy = accuracy_score(y, y_pred)
# accuracy = r2_score(y, y_pred)
print("Accuracy:", round(accuracy, 2))
# print("R2:", round(accuracy, 2))
# print("MAE:", round(np.mean(abs(y - y_pred)), 2))
# plt.scatter(y, y_pred)
# plt.plot(np.linspace(0, 1, 100), np.linspace(0, 1, 100))
# plt.xlim([0, 1])
# plt.ylim([0, 1])

conf_matrix = confusion_matrix(y, y_pred)
print("Confusion Matrix:\n", conf_matrix)
fig, ax = plt.subplots()
sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', cbar=False, ax=ax)
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.title('Confusion Matrix')
# plt.savefig(f'../plots/high_throughput/descriptor_confmat_{system}.png', dpi=150)

# x_range = np.linspace(min(h_alloy_list), max(h_alloy_list), 30)
# y_range = np.linspace(min(h_comp_list), max(h_comp_list), 30)
# xx, yy = np.meshgrid(x_range, y_range)
#
# # Calculate z based on the equation of the plane
# epsilon = 1e-6  # A small threshold to avoid division by zero
# if abs(coef[2]) > epsilon:
# 	zz = -(intercept + coef[0] * xx + coef[1] * yy) / coef[2]
# else:
# 	zz = np.full_like(xx, np.array(temperature_list).mean())
#
# # plt.scatter(h_alloy_list, h_comp_list, c = reduction_flag, cmap = 'coolwarm_r')
# fig = go.Figure(data=[go.Scatter3d(
# 	x=h_alloy_list,
# 	y=h_comp_list,
# 	z=temperature_list,
# 	mode='markers',
# 	marker=dict(
# 		size=5,
# 		color=reduction_flag,  # Set color based on the binary labels
# 		colorscale='plasma',  # Choose a colorscale
# 		opacity=1.0
# 	)
# )])
#
# # Add the decision boundary as a plane
# fig.add_trace(go.Surface(
# 	x=x_range,
# 	y=y_range,
# 	z=zz,
# 	opacity=1,
# 	colorscale='Viridis',
# 	showscale=False
# ))
#
# # Set plot title and labels
# fig.update_layout(
# 	title='3D Scatter Plot with Decision Boundary',
# 	scene=dict(
# 		xaxis_title='H_a',
# 		yaxis_title='H_c',
# 		zaxis_title='$\delta T_{melt}$'
# 	)
# )
#
# # Show the plot
# fig.write_image(f'../plots/high_throughput/descriptor_boundary_{system}.png')
