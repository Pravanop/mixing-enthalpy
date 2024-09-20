import os.path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import mpltern

from calcEnthalpy_package.io.dir_handler import DirHandler
from calcEnthalpy_package.math_operations.thermo_calculations import ThermoMaths
from calcEnthalpy_package.phase_diagram.grid_iterators import GridIterator
from calcEnthalpy_package.io.json_handler import JSONHandler


class TernaryVizualization:
	
	def __init__(self,
				 composition,
				 lattice,
				 meta_data,
				 save_flag,
				 contour_flag):
		
		self.composition = composition
		if len(self.composition) != 3:
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
		self.contour_flag = contour_flag
	
	def find_misc_temperatures(self):
		mol_grid, misc_temp = self.grid_iterator.misc_temperature_across_grid(composition=self.composition,
																			  mol_grid_size=self.mol_grid_size,
																			  lattice=self.lattice,
																			  phase_flag=True,
																			  )
		
		return mol_grid, misc_temp
	
	def find_isotherm(self, temperature):
		mol_grid, stables, _ = self.grid_iterator.e_hull_across_grid(
			composition=self.composition,
			mol_grid_size=self.mol_grid_size,
			lattice=self.lattice,
			single_temp_flag=True,
			temp_gradation=temperature
		)
		return mol_grid, stables
	
	def plot_isotherm(self, temperature):
		mol_grid, stables = self.find_isotherm(temperature)
		t, l, r = mol_grid[:, 0], mol_grid[:, 1], mol_grid[:, 2]
		t = t.reshape(-1, 1)
		l = l.reshape(-1, 1)
		r = r.reshape(-1, 1)
		tm = np.array(stables).reshape(-1, 1)
		fig = plt.figure()
		
		ax = fig.add_subplot(projection="ternary")
		ax.grid()
		cax = ax.inset_axes((1.03, 0.1, 0.05, 0.9), transform=ax.transAxes)
		
		if not self.contour_flag:
			pc = ax.scatter(t, l, r, c=tm, cmap="plasma", marker='h', s=60)
		
		else:
			data = np.concatenate([t, l, r, tm], axis=1)
			
			pc = ax.tricontourf(data[:, 0], data[:, 1], data[:, 2], data[:, 3], cmap="plasma")
		colorbar = fig.colorbar(pc, cax=cax)
		colorbar.set_label('$E_{hull}$ (eV/atom)', rotation=270, va='baseline')
		ax.grid(False)
		ax.set_tlabel(f"{self.composition[0]}")
		ax.set_llabel(f"{self.composition[1]}")
		ax.set_rlabel(f"{self.composition[2]}")
		
		if self.save_flag:
			if self.contour_flag:
				updated_folder_path = DirHandler.mkdir_recursrive(
					folders=['ternary_phase_diagram', "isotherms", "contours", f"{'-'.join(sorted(self.composition, reverse=True))}"], folder_path="../plots")
				fig.savefig(fname=f"{updated_folder_path}{temperature}.png",
							dpi=100)
			else:
				updated_folder_path = DirHandler.mkdir_recursrive(
					folders=['ternary_phase_diagram', "isotherms", "scatters",
							 f"{'-'.join(sorted(self.composition, reverse=True))}"], folder_path="../plots")
				fig.savefig(fname=f"{updated_folder_path}{temperature}.png",
							dpi=100)
		return ax, fig
		
		
	def plot_misc_temperatures(self):
		mol_grid, misc_temp = self.find_misc_temperatures()
		for idx, temp in enumerate(misc_temp):
			if temp == -1:
				misc_temp[idx] = 5000
		
		t, l, r = mol_grid[:, 0], mol_grid[:, 1], mol_grid[:, 2]
		t = t.reshape(-1, 1)
		l = l.reshape(-1, 1)
		r = r.reshape(-1, 1)
		tm = np.array(misc_temp).reshape(-1, 1)
		fig = plt.figure()
		
		ax = fig.add_subplot(projection="ternary")
		ax.grid()
		cax = ax.inset_axes((1.03, 0.1, 0.05, 0.9), transform=ax.transAxes)
		
		if not self.contour_flag:
			pc = ax.scatter(t, l, r, c=tm, cmap="coolwarm", marker = 'h', s = 60)
		
		else:
			data = np.concatenate([t, l, r, tm], axis=1)
			
			pc = ax.tricontourf(data[:, 0], data[:, 1], data[:, 2], data[:, 3], cmap="coolwarm")
		colorbar = fig.colorbar(pc, cax=cax)
		colorbar.set_label('$T_{misc}$ K', rotation=270, va='baseline')
		ax.grid(False)
		ax.set_tlabel(f"{self.composition[0]}")
		ax.set_llabel(f"{self.composition[1]}")
		ax.set_rlabel(f"{self.composition[2]}")
		
		if self.save_flag:
			
			if self.contour_flag:
				updated_folder_path = DirHandler.mkdir_recursrive(
					folders=['ternary_phase_diagram', "misc_temp", "contours"], folder_path="../plots")
				fig.savefig(fname=f"{updated_folder_path}{'-'.join(sorted(self.composition, reverse=True))}.png",
						dpi=100)
			else:
				updated_folder_path = DirHandler.mkdir_recursrive(
					folders=['ternary_phase_diagram', "misc_temp", "scatter"], folder_path="../plots")
				fig.savefig(
					fname=f"{updated_folder_path}{'-'.join(sorted(self.composition, reverse=True))}.png",
					dpi=100)
		return ax, fig

