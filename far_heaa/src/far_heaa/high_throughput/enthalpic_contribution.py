import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.io.json_handler import JSONHandler
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
omega_y_list = []
omega_x_list = []
omega_y_list_misc = []
omega_x_list_misc = []
misc_temp = []
misc_temp_misc = []
tm = ThermoMaths()
# temps = []
temp_reduction = []
giant_combination = []
giant_alloying = []


ele_temp = [tm.avg_T_melt(composition=i, mol_ratio = [0.25, 0.25]) for i in element_list]
sns.histplot(ele_temp)
plt.show()

for key, value in file.items():
    comp_ele = key.split('-')
    for add_ele, temp_list in value.items():
        temp_array = np.array(temp_list)
        temp_array[temp_array==-1000.0]=np.nan
        is_immiscible = temp_list[0] <= 0 # is immiscible


        # print([(system + 1 - x[0])/(system*(system+1))]*system + [x[0]/(system+1)])
        avg_melt_temp = np.array([tm.avg_T_melt(composition=comp_ele + [add_ele], mol_ratio=[(system + 1 - i)/(system*(system+1))]*system + [i/(system+1)]) for i in x])
        temp_array_tot = np.array([avg_melt_temp[idx] - temp_array[idx] if not np.isnan(temp_array[idx]) else np.nan for idx in range(len(temp_array))])
        # print(max(temp_array_tot))
        # temp_reduction.append(temp_array[-1] - temp_array[0])
        if not np.isnan(temp_array_tot[0]) and not np.isnan(temp_array_tot[-1]):
            # temp_reduction.append(avg_melt_temp[-1] - avg_melt_temp[0] - (temp_array_tot[-1] - temp_array_tot[0]))
            temp_reduction.append(temp_array_tot[-1] - temp_array_tot[0])
            # temp_reduction.append(avg_melt_temp[-1] - avg_melt_temp[0])
            omega_x = list(MultinaryCombinations.create_multinary(comp_ele, [2]).values())[0]
            #             # omega_x = omega_x_multi * sum([data[i]['BCC'] for i in omega_x])
            omega_x = sum([data[i]['BCC'] for i in omega_x])
            giant_combination.append(omega_x)

            omega_y = sum([data['-'.join(sorted([i, add_ele]))]['BCC'] for i in comp_ele])
            giant_alloying.append(omega_y)
        # print(np.any(temp_array_tot - temp_array_tot[0]))
        # if np.any(temp_array_tot - temp_array_tot[0]):
        #     temp_reduction.append(max(abs(temp_array_tot - temp_array_tot[0])))
        # else:
        #     temp_reduction.append(max(abs(temp_array_tot - temp_array_tot[0])))




c = ["#009988", "#EE7733", "#0077BB", "#CC3311"]
# print(np.unique(np.array(temp_reduction), return_counts=True))
sns.histplot(temp_reduction, kde = True)
plt.show()


pc = plt.scatter(giant_alloying, giant_combination, c= temp_reduction, zorder = 1)
# plt.plot(temp_reduction)
# plt.scatter(omega_y_list, omega_x_list, c = c[1], zorder = 1)
plt.axhline(y = 0, linestyle='--', c = 'black', zorder = 0)
plt.axvline(x = 0, linestyle='--', c = 'black', zorder = 0)
plt.colorbar(pc)
# plt.legend(['Not Made Misicble', 'Made Miscible'])
plt.xlim([-2, 2])
plt.ylim([-2, 2])
plt.xlabel('Alloying enthalpy')
plt.ylabel('Combination enthalpy')
plt.show()


    #     if is_immiscible:
    #         omega_y = 1/(system+1)**2
    #         # omega_y = 1
    #
    #         indices = np.where(temp_array >= 200)[0] # for safety the misc_T is 200 below avg Tm
    #         if len(indices) != 0 : #there is some alloying that can make it miscible
    #             possible_miscible_keys.append(key)
    #             comp = x[indices]
    #             if key not in results_dict:
    #                 results_dict[key] = [[add_ele, comp]]
    #             else:
    #                 results_dict[key].append([add_ele, comp])
    #
    #
    #             omega_y = sum([data['-'.join(sorted([i, add_ele]))]['BCC'] for i in comp_ele])
    #             omega_y_list.append(omega_y)
    #             omega_x_multi = -(1 + 2 * system) / ((system * (system + 1)) ** 2)
    #             omega_x = list(MultinaryCombinations.create_multinary(comp_ele, [2]).values())[0]
    #             # omega_x = omega_x_multi * sum([data[i]['BCC'] for i in omega_x])
    #             omega_x = sum([data[i]['BCC'] for i in omega_x])
    #             omega_x_list.append(omega_x)
    #             misc_temp.append(temp_array[indices[0]])
    #         else:
    #             omega_y = sum([data['-'.join(sorted([i, add_ele]))]['BCC'] for i in comp_ele])
    #             omega_y_list_misc.append(omega_y)
    #             omega_x_multi = -(1 + 2 * system) / ((system * (system + 1)) ** 2)
    #             omega_x = list(MultinaryCombinations.create_multinary(comp_ele, [2]).values())[0]
    #             # omega_x = omega_x_multi *sum([data[i]['BCC'] for i in omega_x])
    #             omega_x = sum([data[i]['BCC'] for i in omega_x])
    #             omega_x_list_misc.append(omega_x)
    #
    #     else:
    #         miscible_keys.append(key)
    #
    # if key not in miscible_keys and key not in possible_miscible_keys:
    #     never_miscible_keys.append(key)

# fig, ax = plt.subplots(3, 1, figsize = (8, 8))
# print(len(omega_y_list_misc))


# conf_mat = np.zeros((2, 2))
# # trial = [[1, 0], [-1, 2]]
# # print(trial[0][0])
# for i, j in zip(omega_x_list, omega_y_list):
#     if i >=0 and j >=0:
#         conf_mat[0][1] += 1
#     if i <0 and j <0:
#         conf_mat[1][0] += 1
#     if i >=0 and j <0:
#         conf_mat[0][0] += 1
#     if i < 0 and j >=0:
#         conf_mat[1][1] += 1
#
# conf_mat_misc = np.zeros((2, 2))
# for i, j in zip(omega_x_list_misc, omega_y_list_misc):
#     if i >=0 and j >=0:
#         conf_mat_misc[0][1] += 1
#     if i <0 and j <0:
#         conf_mat_misc[1][0] += 1
#     if i >=0 and j <0:
#         conf_mat_misc[0][0] += 1
#     if i < 0 and j >=0:
#         conf_mat_misc[1][1] += 1
#
# per_conf_mat = conf_mat*100/len(omega_y_list)
# per_conf_mat_misc = conf_mat_misc*100/len(omega_y_list_misc)
# print(per_conf_mat)
# print(per_conf_mat_misc)
# print(conf_mat*100/(conf_mat_misc + conf_mat))
#
#
# c = ["#009988", "#EE7733", "#0077BB", "#CC3311"]
# plt.scatter(omega_y_list_misc, omega_x_list_misc, c= c[0], zorder = 1)
# plt.scatter(omega_y_list, omega_x_list, c = c[1], zorder = 1)
# plt.axhline(y = 0, linestyle='--', c = 'black', zorder = 0)
# plt.axvline(x = 0, linestyle='--', c = 'black', zorder = 0)
# plt.legend(['Not Made Misicble', 'Made Miscible'])
# plt.xlim([-2, 2])
# plt.ylim([-2, 2])
# plt.xlabel('Alloying enthalpy')
# plt.ylabel('Combination enthalpy')
# plt.show()

# ax[1].bar(range(len(omega_y_list)), omega_y_list)
# ax[1].bar(range(len(omega_y_list)+1, len(omega_y_list)+len(omega_y_list_misc)+1), omega_y_list_misc)
# ax[1].set_title('Addition Term')
# ax[0].bar(range(len(omega_x_list)), omega_x_list)
# ax[0].bar(range(len(omega_x_list)+1, len(omega_y_list)+len(omega_y_list_misc)+1), omega_x_list_misc)
# ax[0].set_title('Composition Term')
# ax[2].bar(range(len(omega_x_list)), np.array(omega_x_list) + np.array(omega_y_list))
# ax[2].bar(range(len(omega_y_list)+1, len(omega_y_list)+len(omega_y_list_misc)+1), np.array(omega_x_list_misc) + np.array(omega_y_list_misc))
# ax[2].set_title('Total')
# ax[0].set_ylim(-0.2, 0.15)
# ax[1].set_ylim(-0.2, 0.15)
# ax[2].set_ylim(-0.2, 0.15)
# ax[0].set_ylabel('eV/atom')
# ax[1].set_ylabel('eV/atom')
# ax[2].set_ylabel('eV/atom')
# plt.subplots_adjust(hspace = 0.4)
# plt.show()