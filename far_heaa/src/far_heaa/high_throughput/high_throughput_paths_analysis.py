import numpy as np
import matplotlib.pyplot as plt

from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.io.json_handler import JSONHandler
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
for key, value in file.items():
    for add_ele, temp_list in value.items():
        temp_array = np.array(temp_list)
        temp_array[temp_array == -1000.0] = np.nan
        is_immiscible = temp_list[0] <= 0 # is immiscible
        if is_immiscible:
            indices = np.where(temp_array >= 200)[0] # for safety the misc_T is 200 below avg Tm
            if len(indices) != 0 : #there is some alloying that can make it miscible
                possible_miscible_keys.append(key)
                comp = x[indices]
                if key not in results_dict:
                    results_dict[key] = [[add_ele, comp]]
                else:
                    results_dict[key].append([add_ele, comp])

        else:
            miscible_keys.append(key)

    if key not in miscible_keys and key not in possible_miscible_keys:
        never_miscible_keys.append(key)

miscible_keys = list(set(miscible_keys))
never_miscible_keys = list(set(never_miscible_keys))
possible_miscible_keys = list(set(possible_miscible_keys))

assert len(miscible_keys) + len(possible_miscible_keys) + len(never_miscible_keys) == len(list(file.keys()))


print("Miscible: ",len(miscible_keys),
      "\nPossible Miscible by alloying:", len(possible_miscible_keys),
      "\nCan not be made miscible:", len(never_miscible_keys),
      "\nTotal number of alloys:", len(list(file.keys())))


y_max = 1000
y_values_eles = np.linspace(-y_max/1.5, y_max/1.5, len(element_list))
y_values_possible = np.linspace(-y_max, y_max, len(possible_miscible_keys))

# fig, ax = plt.subplots()
# ax.scatter([2]*len(y_values_eles), y_values_eles, c = '#DDAA33', s = 150, edgecolor='black', zorder = 1)
# ax.scatter([1]*len(y_values_possible), y_values_possible, c = '#004488',s = 150, edgecolor='black', zorder = 1)
# for idxk, k in enumerate(element_list):
#     ax.text(s=k, y=y_values_eles[idxk] - 3, x=2.05)

possible_miscible_keys = sorted(possible_miscible_keys)
color_map = plt.get_cmap('Dark2')
count = 0
# element_counter = dict(zip(element_list, [0]*len(element_list)))

binary = list(MultinaryCombinations.create_multinary(element_list, [2]).values())[0]
element_counter = dict(zip(binary, [0]*len(binary)))
for idx, i in enumerate(possible_miscible_keys):
    comp_list = results_dict[i]
    # ax.text(s = i, y = y_values_possible[idx] - 3, x = 1 - system*0.1)
    for j in comp_list:

        x1 = 1
        x2 = 2
        y1 = y_values_possible[idx]
        y2 = y_values_eles[element_list.index(j[0])]
        # ax.plot([x1, x2], [y1, y2], c = color_map(idx), zorder = 0, alpha=0.6)
        count += 1
        ele_list = i.split('-')
        binary_list = list(MultinaryCombinations.create_multinary(ele_list, [2]).values())[0]
        for j in binary_list:
            element_counter[j] += 1

element_counter = dict(sorted(element_counter.items(), key=lambda item: item[1]))
element_counter = {key: value for key, value in element_counter.items() if value != 0}
sns.barplot(x = element_counter.keys(), y = element_counter.values())
plt.xticks(rotation=90)
plt.show()
print(element_counter)
# print(count)
# ax.set_xlim(0.6, 2.4)
# ax.set_xticks([])
# ax.set_yticks([])
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.spines['left'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.set_title('Immiscible Alloys that become miscible after alloying')
# plt.show()
# # plt.savefig("../plots/high_throughput/high_throughput_paths.png", dpi = 200)



