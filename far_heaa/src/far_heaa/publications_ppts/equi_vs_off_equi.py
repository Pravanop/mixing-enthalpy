import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns
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

equi = "/Users/mcube/Desktop/Projects/far_heaa/src/far_heaa/output_data/predictions/equi_Mo-Hf"
off_equi = "/Users/mcube/Desktop/Projects/far_heaa/src/far_heaa/output_data/predictions/non_equi_Mo-Hf"
mpl.rcParams['font.size'] = 12
fig, ax = plt.subplots(figsize=(10, 10))
c = ['black', 'red', 'blue']
s = ['d', 's', 'o']
fig, ax = plt.subplots(4, 1, sharey=True, sharex=True)
c = ["#009988", "#EE7733", "#0077BB", "#CC3311"]
legend = ['Binary', 'Ternary', 'Quarternary', 'Quinary']
for idx, i in enumerate([2, 3, 4, 5]):

    df_offequi = pd.read_csv(f"{off_equi}/misc_T_{i}_BCC.csv")
    df_equi = pd.read_csv(f"{equi}/misc_T_{i}_BCC.csv")
    ax[idx].axvline(x=np.nanmedian(df_offequi['Misc_Temp'] - df_equi['Misc_Temp']),
                    linestyle='--',
                    color = 'black',
                    linewidth=2,
                    zorder = 1)
    sns.histplot(df_offequi['Misc_Temp'] - df_equi['Misc_Temp'],
                 ax=ax[idx],
                 color=c[idx],
                 stat='percent',
                 alpha=0.5)

    # Makes the KDE thicker

    ax[idx].set_xticks([0, 250, 750, 1000])  # Sets the x-tick labels to 0 and 750
    ax[idx].set_xlabel('')  # Sets the x-axis label
    ax[idx].set_ylabel('')  # Sets the x-axis label
    ax[idx].text(s = legend[idx], y = 20, x = 800)


fig.supxlabel('$miscT_{off\_equi}- miscT_{equi}$')
fig.supylabel('Percent')
plt.subplots_adjust(hspace = 0, wspace = 0, bottom = 0.15, left = 0.12, right = 0.95)
plt.savefig("../plots/equi_vs_off_equi.png", dpi = 200)


