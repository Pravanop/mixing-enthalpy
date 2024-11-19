from far_heaa.visualization_toolkit.miscellaneous_visualizations import (
    MiscellaneousVisualizations,
)
from far_heaa.io.metadata_handler import MetadataHandler
import matplotlib.pyplot as plt
mH = MetadataHandler()
meta_data = mH.get_metadata
meta_data["flags"]["equi_flag"] = True
meta_data['file_name']['biased'] = 'bokas_old_omegas_processed'
bD = MiscellaneousVisualizations(meta_data=meta_data, save_flag=True)
#
# ax, fig = bD.plot_convex_hull(composition=["Al", "Fe",], temperature=563)
pdp = bD.plot_convex_hull(composition=["Al", "Fe",], temperature=563)
# # fig = bD.heatmap(element_list=['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf'], lattice = 'BCC')
pdp.show()

# import numpy as np
# from matplotlib import pyplot as plt
# import mpltern
# from far_heaa.grids_and_combinations.grid_creation import CompositionGrid
# from far_heaa.io.metadata_handler import MetadataHandler
# from far_heaa.math_operations.thermo_calculations import ThermoMaths
# from far_heaa.visualization_toolkit.visualizations import Visualizations
# import re
#
# tm = ThermoMaths()
#
# composition = ['Cr', 'Ti']
# transition_temperatures = {
#             "Fe": ["BCC", "FCC", 1180],
#             "Ti": ["HCP", "BCC", 1155],
#             "Hf": ["HCP", "BCC", 2016],
#             "Zr": ["HCP", "BCC", 1136],
#             "Mn": ["BCC", "FCC", 1370],
#         }
# mH = MetadataHandler()
# meta_data = mH.get_metadata
# meta_data['flags']['correction'] = False
# viz = Visualizations(lattice = 'BCC', meta_data=meta_data)
# T = 0
# conv_hull = viz.grid_iterator.convex_hull.make_convex_hull(
#             temperature=T, composition=composition, batch_tag=False
#         )
# conv_hull_T = conv_hull
#
# hsym = tm.calc_multinary_multilattice_mix_enthalpy_sym(
#     composition=composition,
#     transition_temperatures=transition_temperatures,
#     data = viz.data,
#     end_member_dict=viz.end_member,
#     correction = meta_data['flags']['correction'],
#     temperature=T
# )
#
# Ssym = tm.create_Ssym(composition=composition)
#
# Gsym = {}
# for i in hsym:
#     Gsym[i] = tm.G_sym(H_sym=hsym[i], S_sym=Ssym, T = T)
#
# hess = {}
# for i in Gsym:
#     hess[i], xsym = tm.hessian(Gsym[i], composition=composition)
#
#
# mol_grid = CompositionGrid.create_mol_grid(len(composition), 40)
#
# unstable = []
# meta_stable = []
# for mol in mol_grid:
#     mol = list(mol[:len(composition)-1])
#     g_num = np.array([i.subs(list(zip(xsym, mol))) for _, i in Gsym.items()])
#     plt.scatter(mol[0], g_num[0], c = '#BB5566', label = 'BCC')
#     plt.scatter(mol[0], g_num[1], c = '#DDAA33', label = 'FCC')
#     plt.scatter(mol[0], g_num[2], c = '#004488', label = 'HCP')
#
# plt.ylabel('$G_{mix}$', fontsize = 12)
# plt.xlabel('x', fontsize = 12)
# # plt.legend()
# plt.show()






