import numpy as np
from matplotlib import pyplot as plt
import mpltern
from far_heaa.grids_and_combinations.grid_creation import CompositionGrid
from far_heaa.io.metadata_handler import MetadataHandler
from far_heaa.math_operations.thermo_calculations import ThermoMaths
from far_heaa.visualization_toolkit.visualizations import Visualizations
import re

tm = ThermoMaths()

composition = ['Cr', 'V', 'Ti']
transition_temperatures = {
            "Fe": ["BCC", "FCC", 1180],
            "Ti": ["HCP", "BCC", 1155],
            "Hf": ["HCP", "BCC", 2016],
            "Zr": ["HCP", "BCC", 1136],
            "Mn": ["BCC", "FCC", 1370],
        }
mH = MetadataHandler()
meta_data = mH.get_metadata
viz = Visualizations(lattice = 'BCC', meta_data=meta_data)
temp_grid = list(np.linspace(0, 1500, 10))
conv_hull = viz.grid_iterator.temp_iterator(composition=composition, temp_grid=temp_grid)

for T in temp_grid:

    conv_hull_T = conv_hull[T]

    hsym = tm.calc_multinary_multilattice_mix_enthalpy_sym(
        composition=composition,
        transition_temperatures=transition_temperatures,
        data = viz.data,
        end_member_dict=viz.end_member,
        correction = meta_data['flags']['correction'],
        temperature=T
    )

    Ssym = tm.create_Ssym(composition=composition)

    Gsym = {}
    for i in hsym:
        Gsym[i] = tm.G_sym(H_sym=hsym[i], S_sym=Ssym, T = T)

    hess = {}
    for i in Gsym:
        hess[i], xsym = tm.hessian(Gsym[i], composition=composition)

    mol_grid = CompositionGrid.create_mol_grid(len(composition), 40)

    unstable = []
    meta_stable = []
    for mol in mol_grid:
        mol = list(mol[:len(composition)-1])
        g_num = np.array([i.subs(list(zip(xsym, mol))) for _, i in Gsym.items()])
        lattices = list(Gsym.keys())
        lattice = lattices[np.argmin(g_num)]
        H_num = hess[lattice].subs(list(zip(xsym, mol)))
        eigen = tm.find_eigenvalue(H_num)
        is_metastable = True
        is_unstable = False
        for i in eigen:
            if i > 0:
                continue
            if i < 0:
                is_unstable = True
                is_metastable = False

        if is_metastable:
            is_stable = False
            # e_hull =
            meta_stable.append(mol)
        if is_unstable:
            unstable.append(mol)

    fig = plt.figure()

    ax = fig.add_subplot(projection="ternary")
    ax.grid()

    for i in unstable:
        ax.scatter(i[0], i[1],1 - i[0] - i[1], c = '#BB5566')

    for i in meta_stable:
        ax.scatter(i[0],  i[1],  1 - i[0] - i[1], c = '#004488')

    plt.show()






