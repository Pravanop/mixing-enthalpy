import os

import numpy as np
import matplotlib.pyplot as plt
import mpltern
from tqdm import tqdm
from calculateEnthalpy.helper_functions.grid_code import create_mol_grid
from calculateEnthalpy.helper_functions.phaseDiagram import phaseDiagram

pD = phaseDiagram()
composition = "Cr-V-Ti"
phase_diag_dict = pD.make_PD_comp_temp(composition=composition)
temps = list(phase_diag_dict.keys())[::2]
grid_size = 30
mol_grid = create_mol_grid(3, grid_size)
mol_grid = np.array(mol_grid)
ele_list = composition.split('-')
folder_path = "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/plots/phase_diagrams"
os.mkdir(f"{folder_path}/{composition}")
t, l, r = mol_grid[:, 0], mol_grid[:, 1], mol_grid[:, 2]
for temp in tqdm(temps, desc="Creating phase diagrams"):
    stables = []
    for idx, mol in enumerate(mol_grid):
        is_stable = pD.check_stability(mol_ratio=mol,
                                   temp=temp,
                                   conv_hull=phase_diag_dict[temp])
        if is_stable is not None:
            if np.isclose(is_stable[1], 0.0, atol=1e-3):
                stable = 0
            else:
                stable = 1
        else:
            stable = 0
        stables.append(stable)


    ax = plt.subplot(projection="ternary")
    mol_grid = np.array(mol_grid)
    ax.tricontourf(t, l, r, stables)
    ax.grid()
    ax.set_tlabel(f"{ele_list[0]}")
    ax.set_llabel(f"{ele_list[1]}")
    ax.set_rlabel(f"{ele_list[2]}")
    plt.savefig(f"{folder_path}/{composition}/{temp}.png")