import numpy as np
from analyse_cases_functions import plot_everything
from far_heaa.io.json_handler import JSONHandler
from far_heaa.math_operations.thermo_calculations import ThermoMaths
import matplotlib.pyplot as plt
system = 3
element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf']
file = JSONHandler.load_json(folder_path = '/Users/pravanomprakash/Library/CloudStorage/Box-Box/HEA_alloys/done_data/', file_name=f'{system}_add_ele_paths_total_10_BCC_wo_correction')
data = JSONHandler.load_json(folder_path='../database', file_name='bokas_omegas_processed')
mol_grid_size = 5
x = np.linspace(0, 1, mol_grid_size)

results_dict = {}
miscible_keys = []
never_miscible_keys = []
possible_miscible_keys = []
tm = ThermoMaths()
cases = dict(zip(['im-im', 'im-m', 'm-m', 'm-im'], [[],[],[],[]]))
throw_away = []
error_count = 0
counter = {}
for key, value in file.items():
    for add_ele, temp_list in value.items():
        temp_array = np.array(temp_list)
        # temp_array[temp_array == -1000.0] = np.nan

        if np.isnan(temp_array[0]) and np.isnan(temp_array[-1]):
            throw_away.append('.'.join([key, add_ele]))
            continue

        # if np.all(np.isnan(temp_array)):
        #     throw_away.append('.'.join([key, add_ele]))
        #     continue

        if temp_array[0] >= 200:
            indices = list(np.where(temp_array>=200)[0])
            # indices = temp_array[-1] >= 200
            if indices != []:
                cases['m-m'].append('.'.join([key, add_ele]))
            else:
                cases['m-im'].append('.'.join([key, add_ele]))

        if temp_array[0] < 200 or np.isnan(temp_array[0]):
            # indices = temp_array[-1] >= 200
            indices = list(np.where(temp_array >= 200)[0])
            if indices != []:
                speed = min(indices)
                counter.update({speed: counter.get(speed, 0) + 1})
                cases['im-m'].append('.'.join([key, add_ele]))
            else:
                cases['im-im'].append('.'.join([key, add_ele]))



print("Total Cases: ", len(list(file.keys()))*(10-system))
print("Miscible to Miscible: ", len(cases['m-m']))
print('Immiscible to Immiscible: ', len(cases['im-im']) , ' Never immisible:', len(throw_away))
print('Miscible to immiscible: ', len(cases['m-im']))
print('Immiscible to miscible: ', len(cases['im-m']))

counter = dict(sorted(counter.items()))
x = np.linspace(0, 1, mol_grid_size)
import seaborn as sns
import matplotlib

matplotlib.rcParams.update({'font.size': 12})
sns.barplot(x = x[1:]/(system+1), y=list(counter.values()), color='#EE7733')
plt.xlabel('Mole fraction of alloying agent', fontsize=12)
plt.ylabel('Number of alloys made miscible', fontsize=12)
plt.savefig(f'../plots/{system}_imiscible_to_miscible.png')

# assert len(list(file.keys()))*(10-system) == len(cases['m-m']) + len(cases['im-im']) + len(cases['m-im']) + len(cases['im-m']) + len(throw_away)


# for j in ['m-m', 'im-im', 'm-im', 'im-m']:
#     plot_everything(file=file,
#                     data=data,
#                     cases=cases,
#                     case = j,
#                     element_list=element_list,
#                     system=system)
