import json
from typing import List, Union

import numpy as np
from tqdm import tqdm

from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.grids_and_combinations.grid_creation import CompositionGrid
from far_heaa.io.metadata_handler import MetadataHandler
from far_heaa.visualization_toolkit.visualizations import Visualizations

def get_difference(alloy, element_list):
    return list(set(element_list).difference(set(list(alloy.split('-')))))


def find_indices(main_list: List[str], subset: List[str]) -> List[Union[int, None]]:
    indices = []
    for value in subset:
        try:
            index = main_list.index(value)
            indices.append(index)
        except ValueError:
            indices.append(None)
    return indices


element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf']
mol_grid_size = 5
lattice = 'BCC'
for system in [2, 3, 4]:
# system = 2
    alloys = list(MultinaryCombinations.create_multinary(element_list=element_list,
                                                         no_comb=[system+1]).values())[0]

    mH = MetadataHandler()
    meta_data = mH.get_metadata
    meta_data['flags']['correction'] = False
    viz = Visualizations(lattice = lattice,
                             meta_data=meta_data)
    x = list(np.linspace(0, 1, mol_grid_size))
    temp_grid = list(np.linspace(0, 3500, 25))
    data_file = {}
    for idx, alloy in enumerate(tqdm(alloys, desc = f"Iterating through Alloys of dimension {system}")):
        data_file[alloy] = {}
        total_composition = alloy.split('-')
        quarternaries = [[x for i,x in enumerate(total_composition) if i!=j] for j in range(len(total_composition))]
        conv_hull = viz.grid_iterator.temp_iterator(
            composition=total_composition, temp_grid=temp_grid
        )
        for composition in quarternaries:
            N = len(composition)
            n = len(total_composition)
            starting_index = find_indices(
                total_composition, composition
            )
            mol_grid = CompositionGrid.create_high_sym_mol_grid(change_idx=starting_index,
                                                                n=n,
                                                                N=N,
                                                                x=x)
            mol_grid, misc_temp = viz.grid_iterator.misc_temperature_across_grid(
                composition=total_composition,
                mol_grid_size=mol_grid,
                lattice=lattice,
                phase_flag=False,
                conv_hull=conv_hull,
                temp_grid=temp_grid,
                is_differential=True
            )
            ele = list(set(total_composition).difference(composition))[0]
            data_file[alloy][ele] = list(misc_temp[::-1])


    with open(f'./{system}_add_ele_paths_total_{len(element_list)}_{lattice}_wo_correction.json', 'w') as f:
        json.dump(data_file, f, indent=4, ensure_ascii=True)





