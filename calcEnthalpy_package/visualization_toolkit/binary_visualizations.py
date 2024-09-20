import os.path

import numpy as np
from matplotlib import pyplot as plt

from calcEnthalpy_package.io.dir_handler import DirHandler
from calcEnthalpy_package.math_operations.thermo_calculations import ThermoMaths
from calcEnthalpy_package.phase_diagram.grid_iterators import GridIterator
from calcEnthalpy_package.io.json_handler import JSONHandler


class binaryVizualization:
	
	def __init__(self,
				 composition,
				 lattice,
				 meta_data,
				 save_flag):
		
		self.composition = composition
		
		if len(self.composition) != 2:
			raise "Only provide binary compositions!"
		self.lattice = lattice
		self.mol_grid_size = 40
		
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
	
	def find_misc_temperatures(self):
		mol_grid, misc_temp = self.grid_iterator.misc_temperature_across_grid(composition=self.composition,
																			  mol_grid_size=self.mol_grid_size,
																			  lattice=self.lattice,
																			  phase_flag=True,
																			  )
		
		return mol_grid[:-1], misc_temp[:-1]
	
	def plot_misc_temperatures(self):
		mol_grid, misc_temp = self.find_misc_temperatures()
		for idx, temp in enumerate(misc_temp):
			if temp == -1:
				misc_temp[idx] = 5000
		fig, ax = plt.subplots(figsize=(3.2, 3.2))
		ax.set_facecolor('#004488')
		ax.plot(mol_grid[:, 0], misc_temp, color='black', linewidth=2, label='_nolegend_')
		ax.fill(mol_grid[:, 0], misc_temp, color='#DDAA33', label='Stable')
		T_m = [0, 0]
		T_m[1] = self.tm.avg_T_melt(self.composition[0], [1])
		T_m[0] = self.tm.avg_T_melt(self.composition[1], [1])
		ax.plot([0, 1], T_m, '--', color='red', linewidth=1.5)
		
		ax.set_xlim(0, 1)
		ax.set_xticks([0, 1], sorted(self.composition, reverse=True))
		ax.set_ylabel('T (K)')
		ax.set_xlabel('x')
		ax.set_yticks(np.linspace(0, 3000, 6))
		ax.set_ylim(0, max(T_m) + 100)
		ax.set_xticks([0, 1], sorted(self.composition, reverse=True))
		# ax.set_title('-'.join(sorted(self.composition, reverse=True)) + " Phase Diagram")
		plt.subplots_adjust(left=0.22, right=0.95, bottom=0.15)
		ax.legend(['Unstable', 'Average T$_{melt}$'], ncols=2, bbox_to_anchor=(1.1, 1.15), frameon=False)

		if self.save_flag:
			updated_folder_path = DirHandler.mkdir_recursrive(
				folders=['binary_phase_diagram'], folder_path="../plots")
			fig.savefig(
				fname=f"{updated_folder_path}{'-'.join(sorted(self.composition, reverse=True))}.png",
				dpi=100)
			
		return ax, fig
