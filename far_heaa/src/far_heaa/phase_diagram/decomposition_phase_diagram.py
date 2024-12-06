import numpy as np
import matplotlib.pyplot as plt
from far_heaa.visualization_toolkit import visualizations
from far_heaa.io.metadata_handler import MetadataHandler
import re
from far_heaa.grids_and_combinations.grid_creation import CompositionGrid
mH = MetadataHandler()
meta_data = mH.meta_data

meta_data['flags']['im_flag'] = False
viz = visualizations.Visualizations(
    meta_data=meta_data,
    lattice = 'BCC'
)


composition = ['Cr', 'W']

temp_grid = list(np.linspace(0, 2500, 40))
conv_hull = viz.grid_iterator.temp_iterator(composition=composition, temp_grid=temp_grid)


mol_grid = CompositionGrid.create_mol_grid(len(composition), 100)
print(len(mol_grid))

print(conv_hull[0].entries)