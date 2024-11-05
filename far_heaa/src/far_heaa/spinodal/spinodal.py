from matplotlib import pyplot as plt

from far_heaa.io.metadata_handler import MetadataHandler
from far_heaa.visualization_toolkit.visualizations import Visualizations

import numpy as np
from far_heaa.math_operations.spinodal_calculations import calculate_hessian, calculate_eigenvalues


mH = MetadataHandler()
meta_data = mH.get_metadata
viz = Visualizations(lattice= "BCC", meta_data = meta_data)
composition = ['Cr', 'W', 'V']


# Initial Setup
temp_grid = list(np.linspace(0, 2500, 10))
conv_hull = viz.grid_iterator.temp_iterator(composition=composition, temp_grid=temp_grid)
print("Got convex hull")
mol_grid_idx = -1
eigen_value_temp = []
for _, c in conv_hull.items():
    all_entries = c.all_entries
    e_hull_entries = [c.get_e_above_hull(i) for i in all_entries]
    mol_grid = []
    gibbs = []
    for i in all_entries:
        if "BCC" in i.name:
            mol_fraction = [i.composition.get_atomic_fraction(j) for j in composition]
            mol_grid.append(mol_fraction)
            gibbs.append(i.energy / i.composition.num_atoms)

    mol_grid = np.array(mol_grid)[:,[1,0]]
    gibbs = np.array(gibbs)

    if mol_grid_idx != -10:
        mol_fractions = mol_grid[mol_grid_idx]
        hessian = calculate_hessian(mol_fractions, gibbs[mol_grid_idx], mol_grid, gibbs, tolerance_factor=1.0)
        eigenvalues = calculate_eigenvalues(hessian)
        is_stable = np.isclose(e_hull_entries[mol_grid_idx], 0, atol=1e-3, rtol=1e-3)
        eigen_value_temp.append(eigenvalues)

    else:
        for idx, i  in enumerate(mol_grid):
            mol_fractions = i

            hessian = calculate_hessian(mol_fractions, gibbs[idx], mol_grid, gibbs, tolerance_factor=1.0)
            eigenvalues = calculate_eigenvalues(hessian)
            is_stable = np.isclose(e_hull_entries[idx],0, atol = 1e-3, rtol = 1e-3)
            if np.any(eigenvalues<0) and not is_stable:
                print(i)

eigen_value_temp = np.array(eigen_value_temp)

for i in range(len(composition)-1):
    plt.plot(list(conv_hull.keys()), eigen_value_temp[:,i], marker = 'o')

plt.axhline(y = 0)
plt.xticks(list(conv_hull.keys()))
plt.show()


