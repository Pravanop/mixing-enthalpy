import numpy as np
from matplotlib import pyplot as plt
import mpltern
from matplotlib.colors import Normalize
from typing import Tuple

from calcEnthalpy_package.io.dir_handler import DirHandler
from calcEnthalpy_package.math_operations.thermo_calculations import ThermoMaths
from calcEnthalpy_package.phase_diagram.grid_iterators import GridIterator
from calcEnthalpy_package.io.json_handler import JSONHandler


class TernaryVizualization:
	"""
	A class for visualizing ternary phase diagrams and miscibility temperatures for a given alloy composition.

	Args:
		composition (list[str]): A list of three elements in the ternary alloy system.
		lattice (str): The lattice structure (e.g., 'FCC', 'BCC').
		meta_data (dict): A dictionary containing metadata such as grid size, file paths, and flags.
		save_flag (bool): A flag indicating whether to save the visualizations.
		contour_flag (bool): A flag indicating whether to use contour plots or scatter plots.

	Attributes:
		cmap (str): The colormap to use for plotting.
		norm (Normalize): Normalization for color mapping.
		tm (ThermoMaths): An instance of the ThermoMaths class for thermodynamic calculations.
		grid_iterator (GridIterator): An instance of the GridIterator class for grid-based iteration.
		save_flag (bool): A flag indicating whether to save the generated plots.
		contour_flag (bool): A flag indicating whether to generate contour plots.
	"""
	
	def __init__(self,
				 composition: list[str],
				 lattice: str,
				 meta_data: dict,
				 save_flag: bool,
				 contour_flag: bool):
		"""
		Initializes the TernaryVizualization class with the ternary alloy composition, lattice type, and metadata.

		Args:
			composition (list[str]): A list of three elements in the ternary alloy system.
			lattice (str): The lattice structure (e.g., 'FCC', 'BCC').
			meta_data (dict): A dictionary containing metadata such as grid size, file paths, and flags.
			save_flag (bool): A flag indicating whether to save the visualizations.
			contour_flag (bool): A flag indicating whether to use contour plots or scatter plots.

		Raises:
			ValueError: If the provided composition does not contain exactly three elements.
		"""
		self.cmap = 'plasma'
		self.norm = Normalize(vmin=0, vmax=1, clip=False)
		self.composition = composition
		if len(self.composition) != 3:
			raise ValueError("Only provide ternary compositions!")
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
										  flags=meta_data['flags'])
		
		self.save_flag = save_flag
		self.contour_flag = contour_flag
	
	def find_misc_temperatures(self) -> Tuple[np.ndarray, np.ndarray]:
		"""
		Finds miscibility temperatures for the ternary alloy across a mole fraction grid.

		Returns:
			Tuple[np.ndarray, np.ndarray]: A tuple containing:
				- mol_grid (np.ndarray): The mole fraction grid for the ternary alloy.
				- misc_temp (np.ndarray): The miscibility temperatures corresponding to each composition.
		"""
		mol_grid, misc_temp = self.grid_iterator.misc_temperature_across_grid(
			composition=self.composition,
			mol_grid_size=self.mol_grid_size,
			lattice=self.lattice,
			phase_flag=True
		)
		
		return mol_grid, misc_temp
	
	def find_isotherm(self, temperature: float) -> Tuple[np.ndarray, np.ndarray]:
		"""
		Finds the isotherm (temperature slice) of the energy above the convex hull for the ternary alloy.

		Args:
			temperature (float): The temperature in Kelvin to evaluate the isotherm.

		Returns:
			Tuple[np.ndarray, np.ndarray]: A tuple containing:
				- mol_grid (np.ndarray): The mole fraction grid for the ternary alloy.
				- stables (np.ndarray): The stability of the alloy compositions at the given temperature.
		"""
		mol_grid, stables, _ = self.grid_iterator.e_hull_across_grid(
			composition=self.composition,
			mol_grid_size=self.mol_grid_size,
			lattice=self.lattice,
			single_temp_flag=True,
			temp_gradation=temperature
		)
		return mol_grid, stables
	
	def plot_isotherm(self, temperature: float) -> Tuple[plt.Axes, plt.Figure]:
		"""
		Plots the isotherm (temperature slice) of the energy above the convex hull for the ternary alloy.

		Args:
			temperature (float): The temperature in Kelvin to evaluate the isotherm.

		Returns:
			Tuple[plt.Axes, plt.Figure]: A tuple containing the matplotlib axes and figure objects for the plot.
		"""
		mol_grid, stables = self.find_isotherm(temperature)
		
		for idx, stable in enumerate(stables):
			if np.isclose(stable, 0.0, atol=1e-3):
				stables[idx] = 0
			else:
				stables[idx] = 1
		
		t, l, r = mol_grid[:, 0], mol_grid[:, 1], mol_grid[:, 2]
		fig = plt.figure()
		
		ax = fig.add_subplot(projection="ternary")
		ax.grid()
		
		if not self.contour_flag:
			ax.scatter(t, l, r, c=stables, cmap=self.cmap, marker='h', s=60, norm=self.norm)
		else:
			data = np.concatenate([t.reshape(-1, 1), l.reshape(-1, 1), r.reshape(-1, 1), stables.reshape(-1, 1)],
								  axis=1)
			ax.tricontourf(data[:, 0], data[:, 1], data[:, 2], data[:, 3], cmap=self.cmap, norm=self.norm)
		
		ax.grid(False)
		ax.set_tlabel(f"{self.composition[0]}")
		ax.set_llabel(f"{self.composition[1]}")
		ax.set_rlabel(f"{self.composition[2]}")
		
		if self.save_flag:
			if self.contour_flag:
				updated_folder_path = DirHandler.mkdir_recursive(
					folders=['ternary_phase_diagram', "isotherms", "contours",
							 f"{'-'.join(sorted(self.composition, reverse=True))}"], folder_path="../plots")
				fig.savefig(fname=f"{updated_folder_path}{temperature}.png", dpi=100)
			else:
				updated_folder_path = DirHandler.mkdir_recursive(
					folders=['ternary_phase_diagram', "isotherms", "scatters",
							 f"{'-'.join(sorted(self.composition, reverse=True))}"], folder_path="../plots")
				fig.savefig(fname=f"{updated_folder_path}{temperature}.png", dpi=100)
		
		return ax, fig
	
	def plot_misc_temperatures(self) -> Tuple[plt.Axes, plt.Figure]:
		"""
		Plots the miscibility temperatures for the ternary alloy.

		Returns:
			Tuple[plt.Axes, plt.Figure]: A tuple containing the matplotlib axes and figure objects for the plot.
		"""
		mol_grid, misc_temp = self.find_misc_temperatures()
		
		for idx, temp in enumerate(misc_temp):
			if temp == -1:
				misc_temp[idx] = 5000
		
		t, l, r = mol_grid[:, 0], mol_grid[:, 1], mol_grid[:, 2]
		fig = plt.figure()
		
		ax = fig.add_subplot(projection="ternary")
		ax.grid()
		cax = ax.inset_axes((1.03, 0.1, 0.05, 0.9), transform=ax.transAxes)
		
		if not self.contour_flag:
			ax.scatter(t, l, r, c=misc_temp, cmap=self.cmap, marker='h', s=60)
		else:
			data = np.concatenate([t.reshape(-1, 1), l.reshape(-1, 1), r.reshape(-1, 1), misc_temp.reshape(-1, 1)],
								  axis=1)
			ax.tricontourf(data[:, 0], data[:, 1], data[:, 2], data[:, 3], cmap=self.cmap)
		
		colorbar = fig.colorbar(ax, cax=cax)
		colorbar.set_label('$T_{misc}$ K', rotation=270, va='baseline')
		ax.grid(False)
		ax.set_tlabel(f"{self.composition[0]}")
		ax.set_llabel(f"{self.composition[1]}")
		ax.set_rlabel(f"{self.composition[2]}")
		
		if self.save_flag:
			if self.contour_flag:
				updated_folder_path = DirHandler.mkdir_recursive(
					folders=['ternary_phase_diagram', "misc_temp", "contours"], folder_path="../plots")
				fig.savefig(fname=f"{updated_folder_path}{'-'.join(sorted(self.composition, reverse=True))}.png",
							dpi=100)
			else:
				updated_folder_path = DirHandler.mkdir_recursive(
					folders=['ternary_phase_diagram', "misc_temp", "scatter"], folder_path="../plots")
				fig.savefig(fname=f"{updated_folder_path}{'-'.join(sorted(self.composition, reverse=True))}.png",
							dpi=100)
		
		return ax, fig
