import numpy as np
from matplotlib import pyplot as plt

from analyse_cases_functions import plot_everything
from far_heaa.io.json_handler import JSONHandler
from far_heaa.math_operations.thermo_calculations import ThermoMaths

system = 3
element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf']
file = JSONHandler.load_json(folder_path = './done_data/', file_name=f'{system}_add_ele_paths_total_10_BCC_w_im')
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
        composition = key.split('-')
        if add_ele in composition:
            composition.remove(add_ele)

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
                cases['m-im'].append(['-'.join(composition), add_ele])

        if temp_array[0] < 200 or np.isnan(temp_array[0]):
            indices = temp_array[-1] >= 200
            if indices:
                cases['im-m'].append(['-'.join(composition), add_ele, temp_array[-1]])
            else:
                cases['im-im'].append('.'.join([key, add_ele]))


y_ele = element_list
y_ele_values = np.linspace(15, 15 + 2*len(y_ele), len(y_ele))

y_alloy = np.unique(np.array(cases['im-m'])[:,0])
y_alloy_values = np.linspace(0, 0 + 2*len(y_alloy), len(y_alloy))

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random

# Create a directed graph
from pyvis.network import Network
import random

# Create a pyvis Network
net = Network(notebook=True, width="100%", height="750px", bgcolor="#ffffff", font_color="black")

# Add nodes for elements and alloys
for ele in y_ele:
    net.add_node(ele, label=ele, color='lightgreen', shape = 'circle')

for alloy in y_alloy:
    net.add_node(alloy, label=alloy, color='lightblue', shape = 'circle')

connections = [(i[0], i[1], i[2]) for i in np.array(cases['im-m'])]
# connections = sorted(connections, key=lambda x: x[2], reverse=True)
# Add edges with a fixed length but allow dynamic adjustment
for alloy, element, rank in connections:
    net.add_edge(alloy, element, length=100, value= 0.00000001, color = '#EE7733')  # Set a default length (adjust as needed)

# Configure physics for edge length adjustment only
net.toggle_physics(False)
# net.set_options("""
# {
#   "physics": {
#     "enabled": true,
#     "barnesHut": {
#       "centralGravity": 0.0
#     },
#     "stabilization": {
#       "enabled": false
#     }
#   },
#   "edges": {
#     "smooth": false
#   }
# }
# """)

# Display the interactive network
net.show("interactive_graph.html")

