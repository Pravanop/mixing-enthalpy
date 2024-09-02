import json
from typing import Union
from mp_api.client import MPRester
from emmet.core.thermo import ThermoType
import numpy as np
import pandas as pd
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core import Composition
from tqdm import tqdm
import seaborn as sns
from calculateEnthalpy.helper_functions.grid_code import create_multinary, create_mol_grid
from calculateEnthalpy.helper_functions.thermo_math import thermoMaths
import warnings

warnings.filterwarnings("ignore")


class phaseDiagram:
	"""
	Creates a PyMatGen Phase Diagram object for a given phase-space.
	Attributes
	file_path : the path to the json file created from calculateEnthalpyDataset
	binary_dict: the dictionary containing the binary mixing enthalpy values
	tm: Thermodynamic formulae model for calculating properties like enthalpy and entropy
	temp_grid: The temperature grid to calculate the phase diagram for
	"""

	def __init__(self,
				 processed_binary_file_path: str,
				 end_member_file_path: str,
				 grid_size: int,
				 correction: bool,
				 im_flag: bool,
				 equi_flag: bool
				 ) -> None:
		"""

		Args:
			source: a unique source identifier. The naming conventions for the input data can be read in the read me file for this folder.
			lattice: lattice currently working on. Accepted values are 'bcc', 'fcc', 'hcp'
			version_no: the calculateEnthalpyDataset stores folder with a suffix number, that can be used to keep multiple versions
			im_flag: whether the intermetallics have been obtained too
			abs_file_path: path where the data is stored
			binary_dict_path: The folder path for the data folder. Recommended to create a Data/input data folder.
		"""

		self.temp_grid = None
		with open(processed_binary_file_path, "r") as f:
			self.data = json.load(f)

		with open(end_member_file_path, "r") as f:
			self.end_member = json.load(f)

		# self.binary_dict = DataUtils.load_json(folder_path=f"{binary_dict_path}", lattice=lattice, source=source)
		self.tm = thermoMaths()

		self.ele_list = None
		self.grid_size = grid_size

		self.correction = correction
		self.im_flag = im_flag
		self.equi_flag = equi_flag

	def heatmap(self, element_list):
		binaries = list(create_multinary(element_list, no_comb=[2])[2])
		final_data = np.empty((len(element_list), len(element_list)))
		final_data[:] = np.nan
		for pair in binaries:
			print(pair)
			pair_list = pair.split('-')
			if self.correction:
				value = int(self.data[pair]['BCC'] * 1000/4)
			else:
				value = int(self.data[pair]['BCC'] * 1000)
			final_data[element_list.index(pair_list[0])][element_list.index(pair_list[1])] = value
			final_data[element_list.index(pair_list[1])][element_list.index(pair_list[0])] = value

		alloy = pd.DataFrame(final_data, index=element_list, columns=element_list)
		# sns.set(rc={'figure.figsize': (12, 8)})
		fig = sns.heatmap(alloy, cmap="coolwarm", annot=True, fmt='g', mask=alloy.isnull(), square=True, cbar_kws={
			'label': 'Binary Mixing Enthalpy (meV/atom)'}, linewidths=2, linecolor='white')
		return fig

	def upper_limit(self, composition, mol_ratio) -> None:
		avg_tm = self.tm.avg_T_melt(composition, mol_ratio)
		self.temp_grid = np.arange(0, avg_tm + 200, 200, dtype=np.float64)

	def make_convex_hull(self,
						 composition: list,
						 temperature: float):
		n_alloy = len(composition)
		all_combs = create_multinary(element_list=composition, no_comb=list(range(2, n_alloy + 1)))
		pd_entry_input = {}
		pd_entries_list = []
		for dimensionality, alloy_list in all_combs.items():
			if self.im_flag:
				with MPRester("u1TjwfwfTnpF8IolXF9PBY9RT9YauL84", mute_progress_bars=True) as mpr:
					gga = mpr.materials.thermo.search(chemsys=alloy_list,
													  fields=['composition', 'formation_energy_per_atom'],
													  thermo_types=[ThermoType.GGA_GGA_U])
					for i in gga:
						energy = i.formation_energy_per_atom
						name = Composition(i.composition)
						pd_entries_list.append(PDEntry(composition=name,
													   energy=energy * name.num_atoms,
													   name=f'{name.alphabetical_formula}_MP'))

			if self.equi_flag:
				mol_grid = [[1 / dimensionality] * dimensionality]
			else:
				mol_grid = create_mol_grid(int(dimensionality), self.grid_size)
			for alloy_idx, alloy in enumerate(alloy_list):
				# print(type(alloy))
				alloy_list = alloy.split('-')

				for mol_idx, mol_frac in enumerate(mol_grid):
					# print(mol_idx)
					mol_ratio = dict(zip(alloy_list, mol_frac))
					mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}
					mix_enthalpy = self.tm.calc_mutinary_multilattice_mix_Enthalpy(
						mol_ratio=mol_ratio,
						binary_dict=self.data,
						end_member_dict=self.end_member,
						correction=self.correction,
					)

					config_entropy = self.tm.calc_configEntropy(mol_ratio)
					if isinstance(mix_enthalpy, dict):
						for key, enthalpy in mix_enthalpy.items():
							name = Composition(Composition(mol_ratio).get_integer_formula_and_factor()[0])
							pd_entries_list.append(PDEntry(composition=name,
														   energy=self.tm.calc_gibbs_energy(
															   enthalpy,
															   config_entropy,
															   temperature
														   ) * name.num_atoms,
														   name=f'{name.alphabetical_formula}_{key}'))

				for ele in composition:
					for key, value in self.end_member[ele].items():
						name = Composition(ele)

						pd_entries_list.append(PDEntry(composition=name, energy=value * name.num_atoms,
													   name=f'{name.alphabetical_formula}_{key}'))

				for pd_key, value in pd_entry_input.items():
					pd_entries_list.append(PDEntry(composition=pd_key, energy=value))

		phase_diagram = PhaseDiagram(pd_entries_list)
		return phase_diagram

	def make_PD_comp_temp(self,
						  composition: list,
						  temp_grid: list) -> dict[float:PhaseDiagram]:
		"""
		Creates a Phase Diagram object for a temperature grid and composition.
		Args:
			composition: The composition like so - "El1-El2-EL3"

		Returns: A dictionary containing the phase diagram for each temperature for that composition

		"""
		PD_temp_comp_dict = {}

		for idx, temp in enumerate(tqdm(temp_grid, desc="Running Temperature")):
			PD_temp_comp_dict[temp] = self.make_convex_hull(temperature=temp,
															composition=composition)

		return PD_temp_comp_dict

	def check_stability(self,
						mol_ratio: dict,
						temp: float,
						conv_hull: PhaseDiagram,
						mix_enthalpy: float,
						entropy: float) -> Union[tuple[dict, float], None]:
		"""
		Checks whether a composition at that temperature is stable by calculating the energy above hull from the phase diagram
		Args:
			mol_ratio: list of molar ratios for the composition
			temp: Temperature in K
			conv_hull: Phase Diagram object for the phase space

		Returns: A tuple containing decomposition products and energy above hull OR None if composition is deemed stable, or some other error arises.

		"""

		gibbs = self.tm.calc_gibbs_energy(enthalpy=mix_enthalpy,
										  entropy=entropy,
										  temperature=temp)
		pdEntry = self._make_PD_entry(energy=gibbs, mol_ratio=mol_ratio)
		try:

			answer = conv_hull.get_decomp_and_e_above_hull(entry=pdEntry)
			return answer

		except ValueError as e:
			return None

	def find_decomp_products(self,
							 mol_ratio: list,
							 composition: list[str],
							 temperature: float,
							 lattice: str
							 ) -> Union[tuple, int]:
		mol_ratio = dict(zip(composition, mol_ratio))
		mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}
		mix_enthalpy = self.tm.calc_mutinary_multilattice_mix_Enthalpy(mol_ratio=mol_ratio,
																	   binary_dict=self.data,
																	   end_member_dict=self.end_member,
																	   correction=self.correction)[lattice]
		entropy = self.tm.calc_configEntropy(mol_ratio)

		conv_hull = self.make_convex_hull(temperature=float(temperature),
										  composition=composition,
										  )

		is_stable = self.check_stability(mol_ratio=mol_ratio,
										 temp=float(temperature),
										 conv_hull=conv_hull,
										 mix_enthalpy=mix_enthalpy,
										 entropy=entropy)

		if is_stable:
			col1, col2 = [], []
			for key, value in is_stable[0].items():
				col1.append(key.name)
				col2.append(round(value, 3))

			df = pd.DataFrame({'Product': col1, 'Fraction': col2})
			return df, np.round(is_stable[1], 4)
		else:
			return None, 0

	def find_enthalpy_entropy_composition(self,
										  composition,
										  lattice,
										  mol_ratio: list,
										  ):
		mol_ratio = dict(zip(composition, mol_ratio))
		mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}

		if len(mol_ratio.keys()) == 1:
			mix_enthalpy = 0
		else:
			mix_enthalpy = self.tm.calc_mutinary_multilattice_mix_Enthalpy(mol_ratio=mol_ratio,
																		   binary_dict=self.data,
																		   end_member_dict=self.end_member,
																		   correction=self.correction)[lattice]
		entropy = self.tm.calc_configEntropy(mol_ratio)
		return mix_enthalpy, entropy, mol_ratio

	def find_misc_temperature(self,
							  mol_ratio: list,
							  composition: list[str],
							  lattice: str) -> Union[float, None]:
		"""
		Finds miscibility temperature for a given composition
		Args:
			mol_ratio: list of molar ratios for the composition
			composition: The composition like so - "El1-El2-EL3"

		Returns: temperature in K

		"""
		self.upper_limit(composition, mol_ratio)
		mol_ratio = dict(zip(composition, mol_ratio))
		mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}
		mix_enthalpy = self.tm.calc_mutinary_multilattice_mix_Enthalpy(mol_ratio=mol_ratio,
																	   binary_dict=self.data,
																	   end_member_dict=self.end_member,
																	   correction=self.correction)
		if isinstance(mix_enthalpy, dict):
			mix_enthalpy = mix_enthalpy[lattice]

		entropy = self.tm.calc_configEntropy(mol_ratio)
		for idx, temperature in enumerate(self.temp_grid):
			conv_hull = self.make_convex_hull(temperature=float(temperature),
											  composition=composition)

			is_stable = self.check_stability(mol_ratio=mol_ratio,
											 temp=float(temperature),
											 conv_hull=conv_hull,
											 mix_enthalpy=mix_enthalpy,
											 entropy=entropy)

			if is_stable is not None:
				# Using a tolerance of 1 meV/atom for computing energy above the hull
				if np.isclose(is_stable[1], 0.0, atol=1e-3):
					return float(temperature)
				else:
					continue
			else:
				return float(temperature)

		return f"Predicted to melt at {int(self.temp_grid[-1])}"

	@staticmethod
	def _make_PD_entry(mol_ratio: dict,
					   energy: float) -> PDEntry:
		"""
		Helped function to make a PD entry that can be fed into a Phase Diagram.
		Args:
			mol_ratio: list of molar ratios for the composition
			energy: Energy of the phase in eV/atom

		Returns: A PDEntry object

		"""

		name = Composition(Composition(mol_ratio).get_integer_formula_and_factor()[0])
		return PDEntry(composition=name, energy=energy * name.num_atoms)