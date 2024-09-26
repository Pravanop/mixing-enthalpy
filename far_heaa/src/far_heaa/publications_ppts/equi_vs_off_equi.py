import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
from far_heaa.high_throughput.equi_miscT_predictions import EquiMiscTPredictions
from far_heaa.io.metadata_handler import MetadataHandler

# element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Hf', 'Mn', 'Ni']
#
# mH = MetadataHandler()
# meta_data = mH.get_metadata
#
# meta_data['flags']['equi_flag'] = False
# #binaries
# print(meta_data)
# for dim in [2, 3]:
#     equi_miscT = EquiMiscTPredictions(dim = dim,
#                                       element_list = element_list,
#                                       lattice='min',
#                                       meta_data=meta_data,
#                                       save_flag = True,
#                                       generate_plot=False)
#
#     df = equi_miscT.make_predictions()

equi = "/Users/mcube/Desktop/Projects/far_heaa/src/far_heaa/output_data/predictions_equi"
off_equi = "/Users/mcube/Desktop/Projects/far_heaa/src/far_heaa/output_data/predictions"
mpl.rcParams['font.size'] = 16
fig, ax = plt.subplots(figsize=(8, 8))
c = ['black', 'red', 'blue']
s = ['d', 's', 'o']
for idx, i in enumerate([2, 3, 4]):

    df_offequi = pd.read_csv(f"{off_equi}/misc_T_{i}_min.csv")
    df_equi = pd.read_csv(f"{equi}/misc_T_{i}_min.csv")

    ax.scatter(df_equi['Misc_Temp'], df_offequi['Misc_Temp'], c = c[idx], marker = s[idx], alpha = 1.0, s= 60, zorder = 1)

ax.plot(np.linspace(0, 3000), np.linspace(0, 3000), c = 'black', linestyle = '--', zorder = 0)
plt.legend(['Binary', 'Ternary', 'Quaternary'])
plt.xlabel('Equimolar Grid')
plt.ylabel('Off Equimolar Grid')
plt.title('Miscible Temperature Predictions using different grids')
plt.show()
fig.savefig('../plots/equi_vs_off_equi.png', dpi = 100)

