import itertools
import os.path
from typing import Union

import numpy as np
from matplotlib import pyplot as plt

from calcEnthalpy_package.grids_and_combinations import grid_creation
from calcEnthalpy_package.io.dir_handler import DirHandler
from calcEnthalpy_package.io.json_handler import JSONHandler
from calcEnthalpy_package.math_operations.thermo_calculations import ThermoMaths
from calcEnthalpy_package.phase_diagram.grid_iterators import GridIterator


class PathwayVisualizations:
	
	def __init__(self,
				 composition: list,
				 lattice: str,
				 meta_data: dict,
				 save_flag: bool):
		
		self.lattice = lattice
		self.composition = composition
		
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
		t_max = max([self.tm.avg_T_melt(i, mol_ratio=[]) for i in self.composition])
		self.temp_grid = list(np.linspace(0, t_max, 30))
		self.conv_hull = self.grid_iterator.temp_iterator(composition=self.composition,
														  temp_grid=self.temp_grid)
		
		self.mol_grid_size = 5
		self.save_flag = save_flag
	
	def get_rP(self):
		all_pathways = list(itertools.permutations(self.composition, len(self.composition)))
		
		path_dict = {}
		self.x = list(np.linspace(0, 1, self.mol_grid_size))
		for path in all_pathways:
			path = list(path)
			count = 0
			temp_path = []
			
			while count + 1 < len(path):
				start = path[:count + 1]
				end = [path[count + 1]]
				total = start + end
				mol_grid = grid_creation.CompositionGrid.create_high_sym_mol_grid(
					change_idx=list(range(len(total) - 1)),
					x=self.x,
					n=len(total),
					N=len(start)
				)
				mol_grid = mol_grid[::-1]
				_, misc_temp = self.grid_iterator.misc_temperature_across_grid(
					composition=total,
					mol_grid_size=mol_grid,
					lattice=self.lattice,
					phase_flag=True,
					conv_hull=self.conv_hull,
					temp_grid=self.temp_grid,
				)
				for idx, temp in enumerate(misc_temp):
					if temp == -1:
						misc_temp[idx] = 5000
				
				if count == 0:
					misc_temp[0] = self.tm.avg_T_melt(composition=path[:count + 1][0],
													  mol_ratio=[])
				temp_path.append(misc_temp)
				count += 1
			path_dict['-'.join(path)] = temp_path
		return path_dict
	
	@staticmethod
	def get_n_colors_from_cmap(cmap_name, N):
		cmap = plt.get_cmap(cmap_name)  # Get the colormap
		colors = [cmap(i) for i in np.linspace(0, 1, N)]  # Get N evenly spaced colors
		return colors
	
	def text_segregators(self, texts):
		x_dict = {}
		for text in texts:
			if text[1] not in x_dict:
				x_dict[text[1]] = [text]
			else:
				x_dict[text[1]].append(text)
		
		text_dict = {}
		for key, value in x_dict.items():
			value = np.array(value)
			y = value[:, -1].astype(float)
			list_together = self.find_indices_in_range(y, threshold=200)
			print(list_together)
			for i in list_together:
				temp_text = ""
				for j in i:
					temp_text += f"{value[j, 0]}"
					temp_text += ", "
				temp_text = temp_text[:-2]
				text_dict[temp_text] = [float(value[j, 1]) - 0.02*len(temp_text), float(value[j, 2])]
		
		return text_dict
	
	@staticmethod
	def find_indices_in_range(float_list, threshold=200):
		indices_groups = []
		visited = set()  # To track already processed elements
		
		for i in range(len(float_list)):
			if i in visited:
				continue  # Skip elements already part of a group
			
			current_group = [i]  # Start a new group with the current index
			
			# Compare the current element with all subsequent elements
			for j in range(i + 1, len(float_list)):
				if abs(float_list[i] - float_list[j]) <= threshold:
					current_group.append(j)
					visited.add(j)  # Mark the index as visited
			
			# if len(current_group) > 1:
			indices_groups.append(current_group)  # Only keep groups with more than 1 element
		
		return indices_groups
	
	def plot_rP(self):
		
		path_dict = self.get_rP()
		# print(path_dict)
		cmap = self.get_n_colors_from_cmap(cmap_name='coolwarm', N=len(self.composition))
		cmap_dict = dict(zip(self.composition, cmap))
		texts = []
		fig, ax = plt.subplots()
		for key, value in path_dict.items():
			path = key.split('-')
			for idx, i in enumerate(value):
				
				x = np.array(self.x) + idx
				y = np.array(i)
				# ax.scatter(x[0], y[0], s=20, c='black', zorder=1)
				ax.axvline(x=x[0], color='#BBBBBB', alpha=0.3, linestyle='--', zorder=0)
				temp = key.split('-')[:idx + 1]
				texts.append([''.join(sorted(temp)), x[0], y[0]])
				
				if idx == len(value) - 1:
					# ax.scatter(x[-1], y[-1], s=20, c='black', zorder=1)
					texts.append([''.join(sorted(self.composition)), x[-1], y[-1]])
				
				ax.plot(x, y, zorder=0, color=cmap_dict[path[0]])
		
		texts = np.array(texts)
		texts = np.unique(texts, axis=0)
		text_dict = self.text_segregators(texts)
		
		for key, value in text_dict.items():
			t = ax.text(s=key, x=value[0], y=value[1], zorder=100)
			t.set_bbox(dict(facecolor='white', alpha=1.0, pad=2.25, edgecolor='black', linewidth=0.5))
		
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.set_ylim(-100, 4000)
		ax.set_xticks([])
		ax.set_ylabel('Temperature (K)', fontsize=12)
		ax.set_xlabel('Reaction Coordinate', fontsize=12)
		
		if self.save_flag:
			updated_folder_path = DirHandler.mkdir_recursrive(
				folders=['pathways_plots'], folder_path="../plots")
			fig.savefig(f"{updated_folder_path}{'-'.join(sorted(self.composition))}.png", dpi=100)
		
		return ax, fig
