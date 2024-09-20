import numpy as np
import pandas as pd

from calcEnthalpy_package.grids_and_combinations.combination_generation import MultinaryCombinations


class ElementalVisualization:
	def heatmap(self, element_list, genre, sns=None):
		binaries = list(MultinaryCombinations.create_multinary(element_list, no_comb=[2])[2])
		final_data = np.empty((len(element_list), len(element_list)))
		final_data[:] = np.nan
		for pair in binaries:
			pair_list = pair.split('-')
			if self.correction:
				value = int(self.data[pair][genre] * 1000 / 4)
			else:
				value = int(self.data[pair][genre] * 1000)
			final_data[element_list.index(pair_list[0])][element_list.index(pair_list[1])] = value
			final_data[element_list.index(pair_list[1])][element_list.index(pair_list[0])] = value
		
		alloy = pd.DataFrame(final_data, index=element_list, columns=element_list)
		fig = sns.heatmap(alloy, cmap="coolwarm", annot=True, fmt='g', mask=alloy.isnull(), square=True, cbar_kws={
			'label': 'Binary H$_{mix}$ (meV/atom)'}, linewidths=2, linecolor='white')
		return fig