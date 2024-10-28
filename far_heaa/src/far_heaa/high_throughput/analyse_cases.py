import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
from analyse_cases_functions import plot_everything
from calculateEnthalpy.helper_functions.grid_code import create_multinary
from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.io.json_handler import JSONHandler
import seaborn as sns

from far_heaa.math_operations.thermo_calculations import ThermoMaths

system = 3
element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf']
file = JSONHandler.load_json(folder_path = './', file_name=f'{system}_add_ele_paths_total_10_HCP_wo_im')
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
for key, value in file.items():
    for add_ele, temp_list in value.items():
        temp_array = np.array(temp_list)
        temp_array[temp_array == -1000.0] = np.nan

        if np.isnan(temp_array[0]) and np.isnan(temp_array[-1]):
            throw_away.append('.'.join([key, add_ele]))
            continue

        # if np.all(np.isnan(temp_array)):
        #     throw_away.append('.'.join([key, add_ele]))
        #     continue

        if temp_array[0] >= 200:
            # indices = np.where(temp_array>=200)[0]
            indices = temp_array[-1] >= 200
            if indices:
                cases['m-m'].append('.'.join([key, add_ele]))
            else:
                cases['m-im'].append('.'.join([key, add_ele]))

        if temp_array[0] < 200 or np.isnan(temp_array[0]):
            indices = temp_array[-1] >= 200
            if indices:
                cases['im-m'].append('.'.join([key, add_ele]))
            else:
                cases['im-im'].append('.'.join([key, add_ele]))



print("Total Cases: ", len(list(file.keys()))*7)
print("Miscible to Miscible: ", len(cases['m-m']))
print('Immiscible to Immiscible: ', len(cases['im-im']) , ' Never immisible:', len(throw_away))
print('Miscible to immiscible: ', len(cases['m-im']))
print('Immiscible to miscible: ', len(cases['im-m']))

assert len(list(file.keys()))*(10-system) == len(cases['m-m']) + len(cases['im-im']) + len(cases['m-im']) + len(cases['im-m']) + len(throw_away)


plot_everything(file=file,
                data=data,
                cases=cases,
                case = 'm-im',
                element_list=element_list)