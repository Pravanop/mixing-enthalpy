import numpy as np
import pandas as pd
from calculateEnthalpy.helper_functions.grid_code import create_mol_grid

n = 3
mol_grid = create_mol_grid(grid_size=30, n = n)
print(mol_grid.shape)


df = pd.DataFrame(mol_grid, columns=['1', '2', '3'])
df.to_csv(f'./0_{n}space.csv')
print(df.head())