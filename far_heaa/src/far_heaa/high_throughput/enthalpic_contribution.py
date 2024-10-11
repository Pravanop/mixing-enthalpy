import numpy as np
import matplotlib.pyplot as plt

from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.io.json_handler import JSONHandler


system = 3
element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf']
file = JSONHandler.load_json(folder_path = './', file_name=f'{system}_add_ele_paths_total_10')
data = JSONHandler.load_json(folder_path='../database', file_name='bokas_omegas_processed')
mol_grid_size = 5
x = np.linspace(0, 1, mol_grid_size)

results_dict = {}
miscible_keys = []
never_miscible_keys = []
possible_miscible_keys = []
omega_y_list = []
omega_x_list = []
omega_y_list_misc = []
omega_x_list_misc = []
for key, value in file.items():
    comp_ele = key.split('-')
    for add_ele, temp_list in value.items():
        temp_array = np.array(temp_list)
        is_immiscible = temp_list[0] <= 0 # is immiscible

        if is_immiscible:
            omega_y = 1/(system+1)**2
            # omega_y = 1

            indices = np.where(temp_array >= 200)[0] # for safety the misc_T is 200 below avg Tm
            if len(indices) != 0 : #there is some alloying that can make it miscible
                possible_miscible_keys.append(key)
                comp = x[indices]
                if key not in results_dict:
                    results_dict[key] = [[add_ele, comp]]
                else:
                    results_dict[key].append([add_ele, comp])
                omega_y *= sum([data['-'.join(sorted([i, add_ele]))]['BCC'] for i in comp_ele])
                omega_y_list.append(omega_y)
                omega_x_multi = -(1 + 2 * system) / ((system * (system + 1)) ** 2)
                omega_x = list(MultinaryCombinations.create_multinary(comp_ele + [add_ele], [2]).values())[0]
                omega_x = omega_x_multi * sum([data[i]['BCC'] for i in omega_x])
                omega_x_list.append(omega_x)
            else:
                omega_y *= sum([data['-'.join(sorted([i, add_ele]))]['BCC'] for i in comp_ele])
                omega_y_list_misc.append(omega_y)
                omega_x_multi = -(1 + 2 * system) / ((system * (system + 1)) ** 2)
                omega_x = list(MultinaryCombinations.create_multinary(comp_ele + [add_ele], [2]).values())[0]
                omega_x = omega_x_multi * sum([data[i]['BCC'] for i in omega_x])
                omega_x_list_misc.append(omega_x)

        else:
            miscible_keys.append(key)

    if key not in miscible_keys and key not in possible_miscible_keys:
        never_miscible_keys.append(key)

fig, ax = plt.subplots(3, 1, figsize = (8, 8))
print(len(omega_y_list_misc))
ax[1].bar(range(len(omega_y_list)), omega_y_list)
ax[1].bar(range(len(omega_y_list)+1, len(omega_y_list)+len(omega_y_list_misc)+1), omega_y_list_misc)
ax[1].set_title('Addition Term')
ax[0].bar(range(len(omega_x_list)), omega_x_list)
ax[0].bar(range(len(omega_x_list)+1, len(omega_y_list)+len(omega_y_list_misc)+1), omega_x_list_misc)
ax[0].set_title('Composition Term')
ax[2].bar(range(len(omega_x_list)), np.array(omega_x_list) + np.array(omega_y_list))
ax[2].bar(range(len(omega_y_list)+1, len(omega_y_list)+len(omega_y_list_misc)+1), np.array(omega_x_list_misc) + np.array(omega_y_list_misc))
ax[2].set_title('Total')
ax[0].set_ylim(-0.2, 0.15)
ax[1].set_ylim(-0.2, 0.15)
ax[2].set_ylim(-0.2, 0.15)
ax[0].set_ylabel('eV/atom')
ax[1].set_ylabel('eV/atom')
ax[2].set_ylabel('eV/atom')
plt.subplots_adjust(hspace = 0.4)
plt.show()