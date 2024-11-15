import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.io.json_handler import JSONHandler
from far_heaa.math_operations.thermo_calculations import ThermoMaths

# Assuming you have already imported or defined JSONHandler, MultinaryCombinations, and data

element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf']
data = JSONHandler.load_json(folder_path='../database', file_name='bokas_omegas_processed')
mol_grid_size = 5
x = np.linspace(0, 1, mol_grid_size)
lattice = 'BCC'

# Get binaries and corresponding omega values
binaries = list(MultinaryCombinations.create_multinary(element_list=element_list, no_comb=[2]).values())[0]
omegas = [data[i][lattice] for i in binaries]
omegas = sorted(omegas)

# Define colors for each element
element_colors = {element: plt.cm.tab20(i / len(element_list)) for i, element in enumerate(element_list)}

# Plot
fig, ax = plt.subplots(figsize=(10, 8))

# Loop to create each binary bar with two halves (bottom and top colors)
for idx, (binary, omega) in enumerate(zip(binaries, omegas)):
    element1, element2 = binary.split('-')  # Assuming binaries are in the format 'Element1-Element2'
    color1 = element_colors[element1]
    color2 = element_colors[element2]

    # Plot the bottom half of the bar with the color of the first element
    ax.bar(idx, omega / 2, color=color1, width=0.8, edgecolor='black', label=element1 if idx == 0 else "")

    # Plot the top half of the bar with the color of the second element
    ax.bar(idx, omega / 2, bottom=omega / 2, color=color2, width=0.8, edgecolor='black',
           label=element2 if idx == 0 else "")

# Set labels and rotation
ax.set_xticks(range(len(binaries)))
ax.set_xticklabels(binaries, rotation=90)

# Add legend for clarity
handles = [plt.Line2D([0], [0], color=color, lw=4) for color in element_colors.values()]
ax.legend(handles, element_colors.keys(), title="Elements")

plt.show()

# sns.histplot(omegas, bins = len(omegas))
# plt.show()