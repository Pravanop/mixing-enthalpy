import numpy as np
import matplotlib.pyplot as plt
from far_heaa.io.json_handler import JSONHandler
from far_heaa.math_operations.thermo_calculations import ThermoMaths
import seaborn as sns


system = 3
element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf']
file = JSONHandler.load_json(folder_path = './', file_name=f'{system}_add_ele_paths_total_10')
mol_grid_size = 5
x = np.linspace(0, 1, mol_grid_size)

results_dict = {}
miscible_keys = []
never_miscible_keys = []
possible_miscible_keys = []
ele_stats = dict(zip(element_list, [0]*len(element_list)))
tm = ThermoMaths()
temp_diff = []
for key, value in file.items():
    for add_ele, temp_list in value.items():
        temp_array = np.array(temp_list)
        is_immiscible = temp_list[0] <= 0  # is immiscible

        if is_immiscible:
            indices = np.where(temp_array >= 200)[0]
            if len(indices) != 0:
                ele_stats[add_ele] += 1
                temp_diff.append(tm.avg_T_melt(add_ele, [1]) - tm.avg_T_melt(composition=key.split('-'),
                              mol_ratio=[1/system]*system))

sns.histplot(temp_diff, bins=20, color = '#999933', kde = True, zorder = 0)
plt.xlabel('Melting Temp diff between alloying element and alloy (K)')
plt.ylabel('# paths')
plt.show()

ele_stats = dict(sorted(ele_stats.items(), key=lambda item: item[1]))
plt.bar(x = list(ele_stats.keys()), height=list(ele_stats.values()), width = 0.3, color = '#999933', edgecolor = 'black',
        align = 'center')
plt.ylabel('# Alloys made miscible')
plt.xlabel('Alloying element')
plt.xticks(list(ele_stats.keys()))
plt.show()

counter_stats = dict(zip(x, [0]*len(x)))

for key, value in file.items():
    for add_ele, temp_list in value.items():
        temp_array = np.array(temp_list)
        indices = np.where(temp_array >= 200)[0]
        if len(indices) != 0:
            counter_stats[x[min(indices)]] += 1

counter_stats.pop(0.0, None)

plt.bar(x = list(counter_stats.keys()), height=list(counter_stats.values()), width = 0.1, color = '#999933', edgecolor = 'black',
        align = 'center')
plt.ylabel('# Alloys made miscible')
plt.xlabel('Mole Fraction of alloying element')
plt.xticks(list(counter_stats.keys()))
plt.show()
