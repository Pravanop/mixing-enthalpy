import os
import pickle
from statistics import mean
import matplotlib
import numpy as np
from plot_utils import sro_list, array_plotter
import plotly.subplots as sp
from nearest_neighbour import create_neighbor_list
from plotly.subplots import make_subplots
from lookup import enthalpy_model_lookup

composition = 'Cr-W'
model = 'bonds'
nn = 1
atoms = 3456

folder = f"/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/Monte Carlo/dump/{composition}_{atoms}_{model}_{nn}nn/dumps"
# folder = "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/Monte Carlo/dump/Cr-W_2000_enthalpy_1nn/dumps"
lfoldr = os.listdir(folder)
lfoldr.sort(key=lambda x: x.split('_')[1])

temp_range = tuple([x.split('_')[0] for x in lfoldr])

fig_subplots = make_subplots(4, 4,
                             specs=[
                                 [{"type": "scatter3d"}, {"type": "scatter3d"}, {"type": "scatter3d"}, {"type": "scatter3d"}]]*4,
                             subplot_titles=temp_range
                             )
result_fig_traces = []
i, j = 1, 1


lookup, ele_assign = enthalpy_model_lookup(source="pravan",
                                           lattice="bcc",
                                           folder_path="../data/input_data")
for file in lfoldr:

    if file == "plots":
        continue

    if file.endswith(".pickle"):
        with open(f'{folder}/{file}', 'rb') as handle:
            b = pickle.load(handle)

    print(file)
    result = b['structure_trajectory'][-1]
    neighbour_list = create_neighbor_list(result, flag =1)
    sro = sro_list(arr=result, neighbour_list=neighbour_list, ele_assign=ele_assign)
    print(sro)
    result_fig = array_plotter(result, temp = file.split('_')[0], ele_assign = ele_assign)
    print(i,j)
    for trace in range(len(result_fig["data"])):
        fig_subplots.add_trace(result_fig["data"][trace], row = i, col = j)

    if j == 4:
        j = 1
        i += 1
    else:
        j += 1


fig_subplots.show()









