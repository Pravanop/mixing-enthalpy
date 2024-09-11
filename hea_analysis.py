import numpy as np
from matplotlib import pyplot as plt

from calculateEnthalpy.helper_functions.grid_code import create_mol_grid
import matplotlib

matplotlib.rcParams['font.size'] = 16
n = 6
grid_HEA = []
grids = [5, 10,15, 20, 40]
for grid_size in grids:
	print("Grid size: ", grid_size, "Dimensionality: ", n)
	mol_grid = np.round(create_mol_grid(n=n,
										grid_size=grid_size), 4)
	print("Total grid points:", len(mol_grid))

	HEA = 0
	for i in mol_grid:
		if not np.any(i >= 0.50):
			HEA += 1

	grid_HEA.append([len(mol_grid), np.round(HEA/len(mol_grid)*100, 2)])

grid_HEA = np.array(grid_HEA)
plt.plot(grid_HEA[:,0], grid_HEA[:,1], c = 'black', alpha = 0.8, linewidth = 2, marker = 'o', markerfacecolor = 'white', markeredgecolor = 'black', markeredgewidth = 1, ms=10)
plt.title(f"HEAs with order {n}")
plt.ylabel('alloys with no component > 50%')
plt.xlabel('# Grid points')
plt.savefig(f'hea_analysis{n}.png')


