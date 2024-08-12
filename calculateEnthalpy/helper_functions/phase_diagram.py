import json
from typing import Union

import numpy as np
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core import Composition
from tqdm import tqdm

from calculateEnthalpy.helper_functions.data_utils import DataUtils
from calculateEnthalpy.helper_functions.grid_code import create_multinary
from calculateEnthalpy.helper_functions.thermo_math import thermo_maths
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
				 source: str,
				 lattice: str,
				 version_no: int,
				 im_flag: bool,
				 abs_file_path: str,
				 binary_dict_path: str
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
		if not im_flag:
			self.file_path = f"{abs_file_path}/{source}_{lattice}_{version_no}/multinaries_wo_im.json"
		else:
			self.file_path = f"{abs_file_path}/{source}_{lattice}_{version_no}/multinaries_im.json"

		with open(self.file_path, "r") as f:
			self.data = json.load(f)

		self.binary_dict = DataUtils.load_json(folder_path=f"{binary_dict_path}", lattice=lattice, source=source)
		self.tm = thermo_maths(self.binary_dict)
		self.temp_grid = np.arange(0, 3200, 200, dtype=np.float64)
		self.ele_list = None

	def make_PD_composition(self,
							temperature: float,
							composition: str) -> PhaseDiagram:

		"""
		Creates a Phase Diagram object for a given temperature and composition.
		Args:
			temperature: The temperature in K
			composition: The composition like so - "El1-El2-EL3"

		Returns: A pymatgen phase diagram object

		"""

		pd_entries_list = []
		self.ele_list = composition.split('-')
		n_alloy = len(self.ele_list)
		all_combs = create_multinary(element_list=self.ele_list, no_comb=list(range(2, n_alloy + 1)))
		pd_entry_input = {}
		for key, n_nary in all_combs.items():
			#iterating through the various subset dimensionalities present within the alloy starting from 2
			for idx2, subset_comp in enumerate(n_nary):
				#iterating through alloys within that dimension
				temp_subset = self.data[str(key)][subset_comp]
				subset_list = subset_comp.split('-')

				for mol_ratio, value in temp_subset.items():
					#iterating through the mole fraction grid
					if mol_ratio == "intermetallic":
						continue

					#processing mol ratio to get the gibbs energy
					mol_ratio = mol_ratio.split('-')
					mol_ratio = np.array(mol_ratio).astype(float)
					mol_ratio = dict(zip(subset_list, mol_ratio))
					mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}

					name = Composition(Composition(mol_ratio).get_integer_formula_and_factor()[0])
					config_entropy = value['config_entropy']
					pd_entry_input[name] = self.tm.calc_gibbs_energy(
						value['mix_enthalpy'],
						config_entropy,
						temperature
					) * name.num_atoms

				# intermetallics
				if "intermetallic" in temp_subset:
					for idx3, intermetallic in enumerate(temp_subset['intermetallic']):
						name = Composition(intermetallic['formula_pretty'])
						pd_entry_input[name] = intermetallic['formation_energy_per_atom'] * name.num_atoms

				# elements
				for ele in self.ele_list:
					name = Composition(ele)
					pd_entry_input[name] = 0 * name.num_atoms

			pd_entries_list = []
			for pd_key, value in pd_entry_input.items():
				pd_entries_list.append(PDEntry(composition=pd_key, energy=value))

		phase_diagram = PhaseDiagram(pd_entries_list)
		return phase_diagram

	def make_PD_comp_temp(self,
						  composition: str) -> dict[float:PhaseDiagram]:
		"""
		Creates a Phase Diagram object for a temperature grid and composition.
		Args:
			composition: The composition like so - "El1-El2-EL3"

		Returns: A dictionary containing the phase diagram for each temperature for that composition

		"""
		PD_temp_comp_dict = {}
		for idx, temp in enumerate(tqdm(self.temp_grid, desc="Running Temperature")):
			PD_temp_comp_dict[temp] = self.make_PD_composition(temperature=temp,
															   composition=composition)

		return PD_temp_comp_dict

	def check_stability(self,
						mol_ratio: list,
						temp: float,
						conv_hull: PhaseDiagram) -> Union[tuple[dict, float], None]:
		"""
		Checks whether a composition at that temperature is stable by calculating the energy above hull from the phase diagram
		Args:
			mol_ratio: list of molar ratios for the composition
			temp: Temperature in K
			conv_hull: Phase Diagram object for the phase space

		Returns: A tuple containing decomposition products and energy above hull OR None if composition is deemed stable, or some other error arises.

		"""
		mol_ratio = dict(zip(self.ele_list, mol_ratio))
		mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}
		comp = '-'.join(list(mol_ratio.keys()))
		mix_enthalpy = self.tm.calc_multinary_mixEnthalpy(alloy_comp=comp, mol_ratio=mol_ratio)
		entropy = self.tm.calc_configEntropy(mol_ratio)
		gibbs = self.tm.calc_gibbs_energy(enthalpy=mix_enthalpy,
										  entropy=entropy,
										  temperature=temp)
		pdEntry = self._make_PD_entry(energy=gibbs, mol_ratio=mol_ratio)

		try:
			return conv_hull.get_decomp_and_e_above_hull(entry=pdEntry)

		except:
			return None

	def find_misc_temperature(self,
							  mol_ratio: list,
							  composition: str) -> Union[float, None]:
		"""
		Finds miscibility temperature for a given composition
		Args:
			mol_ratio: list of molar ratios for the composition
			composition: The composition like so - "El1-El2-EL3"

		Returns: temperature in K

		"""

		for idx, temperature in enumerate(self.temp_grid):
			conv_hull = self.make_PD_composition(temperature=float(temperature),
												 composition=composition)
			is_stable = self.check_stability(mol_ratio=mol_ratio, temp=float(temperature), conv_hull=conv_hull)
			if is_stable is not None:
				# Using a tolerance of 1 meV/atom for computing energy above the hull
				if np.isclose(is_stable[1], 0.0, atol=1e-3):
					return float(temperature)
				else:
					return None

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
