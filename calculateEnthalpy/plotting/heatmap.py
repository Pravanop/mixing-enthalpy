import matplotlib.pyplot as plt
import numpy as np
from data_utils import load_json
import pandas as pd
import seaborn as sns
sns.set_theme(style = 'dark')
lattice = "bcc"
source = "Cr_W"
binary = load_json(folder_path ="../data/input_data/",
                   source = source,
                   lattice = lattice)
keys = list(binary.keys())
element_list = ['Cr', 'V', 'Ti', 'W', 'Ta', 'Hf', 'Fe', 'Mn']
element_list.sort()
final_data = np.empty((len(element_list), len(element_list)))
final_data[:] = np.nan
for pair in keys:
	pair_list = pair.split('-')
	final_data[element_list.index(pair_list[0])][element_list.index(pair_list[1])] = int(binary[pair]*1000)
	final_data[element_list.index(pair_list[1])][element_list.index(pair_list[0])] = int(binary[pair]*1000)
alloy = pd.DataFrame(final_data, index = element_list, columns = element_list)
sns.heatmap(alloy, cmap = "coolwarm", annot = True,  fmt='g', mask=alloy.isnull(), square = True, cbar_kws={
		'label': 'Binary Mixing Enthalpy (meV/atom)'}, linewidths = 2, linecolor = 'white')
plt.show()


	