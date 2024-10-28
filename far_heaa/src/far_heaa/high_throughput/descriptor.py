import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix
from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.io.json_handler import JSONHandler
import seaborn as sns

from far_heaa.math_operations.thermo_calculations import ThermoMaths

system = 3
element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf']
file = JSONHandler.load_json(folder_path = './', file_name=f'{system}_add_ele_paths_total_10_wo_im')
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
temperature_list = []
reduction = []
for key, value in file.items():
    for add_ele, temp_list in value.items():

        temp_array = np.array(temp_list)
        temp_array[temp_array == -1000.0] = np.nan

        if np.isnan(temp_array[0]) or np.isnan(temp_array[-1]):
            continue

        h_alloy = 0
        for j in key.split('-'):
            h_alloy += data['-'.join(sorted([j, add_ele]))]['BCC']
        h_alloy_list.append(h_alloy)

        h_comp = 0
        for j in list(MultinaryCombinations.create_multinary(key.split('-'), no_comb=[2]).values())[0]:
            h_comp += data[j]['BCC']
        h_comp_list.append(h_comp)

        alloy_temp = tm.avg_T_melt(composition=key.split('-'), mol_ratio=[0.33, 0.33, 0.34])
        element_temp = tm.avg_T_melt(composition=key.split('-') + [add_ele], mol_ratio=[0.25, 0.25, 0.25, 0.25])
        temperature_list.append(float(alloy_temp - element_temp))
        if np.isnan(temp_array[-1]):
            reduction.append(-5000 - temp_array[0])
        elif np.isnan(temp_array[0]):
            reduction.append(temp_array[-1] + 5000)
        else:
            reduction.append(temp_array[-1] - temp_array[0])


reduction_flag = [0 if i < 0 else 1 for i in reduction]

df = pd.DataFrame([h_alloy_list, h_comp_list, temperature_list, reduction_flag]).T
df.columns = ['Ha', 'Hc', 'delta meltT', 'Miscibility']

fig, axs = plt.subplots(3, 1, sharey=True)

# Plotting each feature in a separate subplot
features = ['Ha', 'Hc', 'delta meltT']
for i, feature in enumerate(features):
    sns.histplot(df[feature], bins=15, color='#BB5566', edgecolor='black', stat = 'percent', ax = axs[i])
    axs[i].set_xlabel(feature)
    axs[i].set_ylabel('Frequency')

plt.subplots_adjust(wspace=0 , hspace = 0)
# Adjust layout for better spacing
plt.tight_layout()
plt.show()

pairplot = sns.pairplot(df, hue='Miscibility')
legend = pairplot._legend
legend.set_title("Miscibility")  # Set a custom title
legend.get_texts()[0].set_text("Increases")  # Change the text for the first group
legend.get_texts()[1].set_text("Decreases")
plt.show()
sns.heatmap(df.corr(), annot=True,
            xticklabels=['Ha', 'Hc', 'delta meltT', 'Miscibility'],
            yticklabels=['Ha', 'Hc', 'delta meltT', 'Miscibility'])
print(df.info())
plt.show()

model = LogisticRegression()
X = df[['Ha', 'Hc', 'delta meltT']].to_numpy()
y = df['Miscibility'].to_numpy()
print(X.shape, y.shape)
model.fit(X, y)
coef = model.coef_[0]
intercept = model.intercept_[0]

# Example prediction
y_pred = model.predict(X)
accuracy = accuracy_score(y, y_pred)
print("Accuracy:", round(accuracy, 2))

conf_matrix = confusion_matrix(y, y_pred)
print("Confusion Matrix:\n", conf_matrix)
sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', cbar=False)
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.title('Confusion Matrix')
plt.show()
x_range = np.linspace(min(h_alloy_list), max(h_alloy_list), 30)
y_range = np.linspace(min(h_comp_list), max(h_comp_list), 30)
xx, yy = np.meshgrid(x_range, y_range)

# Calculate z based on the equation of the plane
epsilon = 1e-6  # A small threshold to avoid division by zero
if abs(coef[2]) > epsilon:
    zz = -(intercept + coef[0] * xx + coef[1] * yy) / coef[2]
else:
    zz = np.full_like(xx, np.array(temperature_list).mean())

# plt.scatter(h_alloy_list, h_comp_list, c = reduce_flag, cmap = 'coolwarm_r')
fig = go.Figure(data=[go.Scatter3d(
    x=h_alloy_list,
    y=h_comp_list,
    z=temperature_list,
    mode='markers',
    marker=dict(
        size=5,
        color=reduction_flag,  # Set color based on the binary labels
        colorscale='viridis',  # Choose a colorscale
        opacity=1.0
    )
)])

# Add the decision boundary as a plane
fig.add_trace(go.Surface(
    x=x_range,
    y=y_range,
    z=zz,
    opacity=1,
    colorscale='Viridis',
    showscale=False
))

# Set plot title and labels
fig.update_layout(
    title='3D Scatter Plot with Decision Boundary',
    scene=dict(
        xaxis_title='H_a',
        yaxis_title='H_c',
        zaxis_title='$\delta T_{melt}$'
    )
)

# Show the plot
fig.show()