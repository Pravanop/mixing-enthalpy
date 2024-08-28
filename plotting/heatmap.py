from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
from calculateEnthalpy.helper_functions.data_utils import  DataUtils
import pandas as pd
import seaborn as sns
sns.set_theme(style = 'dark')
lattice = "hcp"
source = "bokas"
binary = DataUtils.load_json(folder_path ="../data/input_data/",
				   source = source,
				   lattice = lattice)
keys = list(binary.keys())
element_list = ['Fe', 'Co', 'Ni', 'Cu','Mn','Mg','Al','Ru', 'Rh', 'Pd', 'Pt', 'Ag','W','Cr', 'Ga', 'In', 'Sn']
element_list.sort()
element_combs = ["-".join(sorted(list(i))) for i in list(combinations(element_list, 2))]
final_data = np.empty((len(element_list), len(element_list)))
final_data[:] = np.nan
# for ele in element_list:
for pair in element_combs:
	pair_list = pair.split('-')
	print(pair_list)
	final_data[element_list.index(pair_list[0])][element_list.index(pair_list[1])] = int(binary[pair]*1000)
	final_data[element_list.index(pair_list[1])][element_list.index(pair_list[0])] = int(binary[pair]*1000)
alloy = pd.DataFrame(final_data, index = element_list, columns = element_list)
sns.set (rc = {'figure.figsize':(12, 8)})
sns.heatmap(alloy, cmap = "coolwarm", annot = True,  fmt='g', mask=alloy.isnull(), square = True, cbar_kws={
		'label': 'Binary Mixing Enthalpy (meV/atom)'}, linewidths = 2, linecolor = 'white')
plt.savefig(f"../plots/heatmap_alloy_Ga_{lattice}.png", dpi=300)
plt.show()


	