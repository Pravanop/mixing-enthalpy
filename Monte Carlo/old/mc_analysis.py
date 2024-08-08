from statistics import mean
import matplotlib
import numpy as np
from plot_utils import retrieve_latest_dump, sro_list, array_plotter
from nearest_neighbour import create_neighbor_list

font = {
    'family': 'Helvetica',
    'size': 16
}
matplotlib.rc('font', **font)

b, path = retrieve_latest_dump(filter='V-W')
energy = b['energy_trajectory'] / 1000
structure = b['structure_trajectory']
split = np.where(energy == 0.1)[0][0]

structure_warm = structure[:split]
structure_main = structure[split + 1:]
# ele_dict = {
#     1: 0.2,
#     2: 0.2,
#     3: 0.2,
#     4: 0.2,
#     5: 0.2
# }
ele_dict = {2:0.5,
            3:0.5}
arr = structure_main[-1]

indices = np.indices(dimensions=arr.shape).reshape(3, -1).T
neighbour_list = create_neighbor_list(arr, flag =2)
path_parts = path.split('_')
# array_plotter(structure_warm[0])
# array_plotter(structure_warm[-1])
array_plotter(structure_warm[0], f'Composition: {path_parts[0]}, No of Atoms: {path_parts[1]}, Temperature: {path_parts[2]}')

sro = sro_list(arr = arr, neighbour_list=neighbour_list, ele_dict=ele_dict)

for key,values in sro.items():
    print(key, np.round(mean(values),3))