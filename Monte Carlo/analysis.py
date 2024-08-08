import os
import pickle
from statistics import mean
import matplotlib
import numpy as np
from plot_utils import retrieve_latest_dump, sro_list, array_plotter
import plotly.subplots as sp
from nearest_neighbour import create_neighbor_list
from plotly.subplots import make_subplots
# dump_dict, path = retrieve_latest_dump(filter='V-W')
composition = 'Cr-W'
model = 'bonds'
nn = 2
atoms = 2000
folder = f"/Users/pravanomprakash/Documents/Projects/highEntropyAlloys/Monte Carlo/dump/{composition}_{atoms}_{model}_{1}nn/dumps"
lfoldr = os.listdir(folder)
lfoldr.sort(key=lambda x: x.split('_')[1])

temp_range = tuple([x.split('_')[0] for x in lfoldr])

fig_subplots = make_subplots(4, 4,
                             specs=[
                                 [{"type": "scatter3d"}, {"type": "scatter3d"}, {"type": "scatter3d"}, {"type": "scatter3d"}],
                                 [{"type": "scatter3d"}, {"type": "scatter3d"}, {"type": "scatter3d"}, {"type": "scatter3d"}],
[{"type": "scatter3d"}, {"type": "scatter3d"}, {"type": "scatter3d"}, {"type": "scatter3d"}],
[{"type": "scatter3d"}, {"type": "scatter3d"}, {"type": "scatter3d"}, {"type": "scatter3d"}]],
                             subplot_titles=temp_range
                             )
result_fig_traces = []
i, j = 1, 1


for file in lfoldr:

    if file == "plots":
        continue

    if file.endswith(".pickle"):
        with open(f'{folder}/{file}', 'rb') as handle:
            b = pickle.load(handle)

    print(file)
    result = b['structure_trajectory'][-1]

    result_fig = array_plotter(result, temp = file.split('_')[0])
    print(i,j)
    for trace in range(len(result_fig["data"])):
        fig_subplots.add_trace(result_fig["data"][trace], row = i, col = j)

    if j == 4:
        j = 1
        i += 1
    else:
        j += 1


fig_subplots.show()







