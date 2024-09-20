import os
from typing import Union
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from calcEnthalpy_package.io.dir_handler import DirHandler
from calcEnthalpy_package.io.json_handler import JSONHandler
from calcEnthalpy_package.math_operations.thermo_calculations import ThermoMaths
from calcEnthalpy_package.phase_diagram.grid_iterators import GridIterator
from calcEnthalpy_package.grids_and_combinations import grid_creation


class MatrixHeatmap:
	
	def __init__(self,
				 composition: list,
				 add_ele: Union[list, str],
				 lattice: str,
				 meta_data: dict,
				 save_flag: bool,
				 type: str = 'add'):
		self.composition = composition
		if isinstance(add_ele, str):
			self.add_ele = [add_ele]
		else:
			self.add_ele = add_ele
		
		self.total_composition = self.composition + self.add_ele
		self.total_composition = list(set(self.total_composition))
		
		if type == 'add':
			self.N = len(self.composition)
			self.n = len(self.total_composition)
			self.n_alloy = len(self.composition)
			self.starting_index = self.find_indices(self.total_composition, self.composition)
			print(self.starting_index)
		elif type == 'transmutate':
			assert len(self.add_ele) == 2
			assert self.add_ele[0] in self.composition
			self.N = len(self.add_ele)
			self.n = len(self.total_composition)
			self.n_alloy = len(self.composition)
			self.transmutation_indice = self.find_indices(self.total_composition, self.add_ele)
			self.end_composition = list(np.copy(self.total_composition))
			self.end_composition.remove(self.add_ele[0])
		
		self.type = type
		self.lattice = lattice
		self.mol_grid_size = 10
		
		self.tm = ThermoMaths()
		grid_size = meta_data['grid_size']
		
		if meta_data['flags']['correction']:
			data = JSONHandler.load_json(folder_path=meta_data['folder_path'],
										 file_name=meta_data['file_name']['biased'])
		else:
			data = JSONHandler.load_json(folder_path=meta_data['folder_path'],
										 file_name=meta_data['file_name']['unbiased'])
		
		end_member = JSONHandler.load_json(folder_path=meta_data['folder_path'], file_name=meta_data['end_member'])
		
		self.grid_iterator = GridIterator(grid_size=grid_size,
										  tm=self.tm,
										  data=data,
										  end_member=end_member,
										  api_key=meta_data['api_key'],
										  flags=meta_data['flags']
										  )
		
		self.save_flag = save_flag
	
	@staticmethod
	def find_indices(main_list, subset):
		indices = []
		for value in subset:
			try:
				index = main_list.index(value)  # Find the index of the value in the main list
				indices.append(index)
			except ValueError:
				indices.append(None)  # If the value is not found, append None or handle as needed
		return indices
	
	def get_ehull_matrix(self):
		
		self.x = list(np.linspace(0, 1, self.mol_grid_size))
		
		if self.type == 'add':
			mol_grid = grid_creation.CompositionGrid.create_high_sym_mol_grid(
				change_idx=self.starting_index,
				x=self.x,
				n=self.n,
				N=self.N
			)
		elif self.type == 'transmutate':
			mol_grid = grid_creation.CompositionGrid.create_mol_grid_transmutation(
				x = self.x,
				n=self.n,
				transmutation_indice=self.transmutation_indice
			)
		# print(mol_grid)
		mol_grid = mol_grid[::-1]
		
		mol_grid, e_hulls, temp_grid = self.grid_iterator.e_hull_across_grid(
			composition=self.total_composition,
			lattice=self.lattice,
			mol_grid_size=mol_grid,
			single_temp_flag=False,
			temp_gradation=400
		)
		return mol_grid, e_hulls, temp_grid
	
	def plot_ehull_matrix(self):
		
		mol_grid, e_hulls, temp_grid = self.get_ehull_matrix()
		
		df = pd.DataFrame(e_hulls)
		df.columns = temp_grid
		
		df = df.T
		df = df.apply(lambda x: x * 1000)
		positions = []
		for col in df.columns:
			try:
				first_zero_index = df[df[col] == 0].index[0]
			except IndexError:
				first_zero_index = df.index[-1]
			
			positions.append(first_zero_index)
		
		cmap = sns.cubehelix_palette(start=.5, rot=-.61, light=.98, dark=.35, hue=1, as_cmap=True)
		sns.set_theme(rc={'figure.figsize': (6.4, 6.4)})
		sns.set(font_scale=1.2)
		fig, ax = plt.subplots()
		g = sns.heatmap(df,
						yticklabels=np.array(temp_grid).astype(int),
						cmap=cmap,
						# xticklabels=np.round(mol_grid[:,-1],2),
						xticklabels=np.round(self.x,2),
						cbar_kws={'label': '$E_{hull}$ (meV/atom)', "shrink": 0.8},
						ax=ax)
		g.set_yticklabels(g.get_yticklabels(), rotation=0)
		g.set_xticklabels(g.get_xticklabels(), rotation=0)
		ax.axes.invert_yaxis()
		ax.axhline(y=0, color='k', linewidth=3)
		ax.axhline(y=df.shape[0], color='k', linewidth=3)
		ax.axvline(x=0, color='k', linewidth=3)
		ax.axvline(x=10, color='k', linewidth=3)
		
		count_prev = np.where(temp_grid == positions[0])[0][0] + 1
		for idx, i in enumerate(positions):
			idx2 = np.where(temp_grid == i)[0][0] + 1
			if idx2 != count_prev:
				ax.plot([idx, idx], [count_prev, idx2], color='black', linestyle='--')
				count_prev = idx2
			ax.plot([idx, idx + 1], [idx2, idx2], color='black', linestyle='--')
		ax.text(s='-'.join(self.composition), y=-1, x=-1)
		if self.type == 'transmutate':
			ax.text(s='-'.join(self.end_composition), y=-1, x=len(mol_grid)-0.3)
		else:
			ax.text(s='-'.join(self.total_composition), y=-1, x=len(mol_grid) - 0.3)
		ax.set_xlabel("X")
		ax.set_ylabel("T (K)")
		plt.subplots_adjust(bottom=0.15, top=0.9, left=0.14, right=0.98)

		if self.save_flag:
			if self.type == 'transmutate':
				updated_folder_path = DirHandler.mkdir_recursrive(
					folders=['heatmap_plots', "transmutate"], folder_path="../plots")
			
				fig.savefig(
					fname=f"{updated_folder_path}{''.join(self.composition)}_{''.join(self.end_composition)}.png",
					dpi=100)
			elif self.type == 'add':
				updated_folder_path = DirHandler.mkdir_recursrive(
					folders=['heatmap_plots', "add"], folder_path="../plots")
				fig.savefig(
					fname=f"{updated_folder_path}{''.join(self.total_composition)}_{''.join(self.composition)}.png",
					dpi=100)
		return ax, fig
