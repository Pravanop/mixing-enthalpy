import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Composition
from calcEnthalpy_package.grids_and_combinations.combination_generation import MultinaryCombinations
from calcEnthalpy_package.grids_and_combinations.grid_creation import CompositionGrid
from calcEnthalpy_package.data_processing.intermetallic_extraction import IntermetallicExtractions
from calcEnthalpy_package.phase_diagram.pdEntry_local import PDEntryLocal
from calcEnthalpy_package.math_operations.thermo_calculations import ThermoMaths


class ConvexHull:
	
	def __init__(self,
				 end_member: dict,
				 data: dict,
				 flags: dict,
				 tm: ThermoMaths,
				 api_key: str,
				 grid_size: int
				 ):
		
		self.transition_temperatures = {
			'Fe': ['BCC', 'FCC', 1180],
			'Ti': ['HCP', 'BCC', 1155],
			'Hf': ['HCP', 'BCC', 2016],
			'Zr': ['HCP', 'BCC', 1136],
			'Mn': ['BCC', 'FCC', 1370]
		}
		self.end_member = end_member
		self.correction = flags['correction']
		self.data = data
		self.tm = tm
		self.im_flag = flags['im_flag']
		self.equi_flag = flags['equi_flag']
		self.api_key = api_key
		self.grid_size = grid_size
	
	def make_convex_hull(self,
						 composition: list,
						 temperature: float,
						 batch_tag: bool = False,
						 **kwargs) -> PhaseDiagram:
		
		pd_entries_list = []
		
		if batch_tag and self.im_flag:
			if kwargs['im']:
				pd_entries_list += kwargs['im']
		
		n_alloy = len(composition)
		all_combs = MultinaryCombinations.create_multinary(element_list=composition,
														   no_comb=list(range(2, n_alloy + 1)))
		
		for dimensionality, alloy_list in all_combs.items():
			if not batch_tag and self.im_flag:
				pd_entries_list += IntermetallicExtractions.get_MP_intermetallic(alloy_list=alloy_list,
																				 api_key=self.api_key)
			if self.equi_flag:
				mol_grid = [[1 / dimensionality] * dimensionality]
			else:
				mol_grid = CompositionGrid.create_mol_grid(int(dimensionality), self.grid_size)
			
			for alloy_idx, alloy in enumerate(alloy_list):
				alloy_list = alloy.split('-')
				
				for mol_idx, mol_frac in enumerate(mol_grid):
					mol_ratio = dict(zip(alloy_list, mol_frac))
					mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}
					
					if len(mol_ratio.keys()) == 1:
						continue
					
					mix_enthalpy = self.tm.calc_mutinary_multilattice_mix_Enthalpy(
						mol_ratio=mol_ratio,
						binary_dict=self.data,
						end_member_dict=self.end_member,
						correction=self.correction,
						temperature=temperature,
						transition_temperatures=self.transition_temperatures,
						model='regular'
					)
					config_entropy = self.tm.calc_configEntropy(mol_ratio)
					
					if isinstance(mix_enthalpy, dict):
						for key, enthalpy in mix_enthalpy.items():
							name = Composition(Composition(mol_ratio).get_integer_formula_and_factor()[0])
							pd_entries_list.append(PDEntryLocal(composition=name,
																energy=self.tm.calc_gibbs_energy(
																	enthalpy,
																	config_entropy,
																	temperature
																) * name.num_atoms,
																name=f'{name.alphabetical_formula}_{key}'))
		
		pd_entries_list = self.process_elements_hull(composition, temperature, pd_entries_list)
		conv_hull = PhaseDiagram(pd_entries_list)
		
		return conv_hull
	
	def update_convex_hull(self, composition, conv_hull, temperature):
		pd_entries = conv_hull.all_entries
		new_pd_entries = []
		for idx, entry in enumerate(pd_entries):
			lattice = entry.name.split('_')[-1]
			if 'MP' in entry.name:
				new_pd_entries.append(entry)
			else:
				mol_ratio = {el.symbol: np.round(entry.composition.get_atomic_fraction(el), 4) for el in
							 entry.composition.elements}
				if len(list(mol_ratio.keys())) > 1:
					mix_enthalpy = self.tm.calc_mutinary_multilattice_mix_Enthalpy(
						mol_ratio=mol_ratio,
						binary_dict=self.data,
						end_member_dict=self.end_member,
						correction=self.correction,
						temperature=temperature,
						transition_temperatures=self.transition_temperatures,
						model='regular'
					)
					
					config_entropy = self.tm.calc_configEntropy(mol_ratio)
					if isinstance(mix_enthalpy, dict):
						entry.update_energy(self.tm.calc_gibbs_energy(
							mix_enthalpy[lattice],
							config_entropy,
							temperature
						) * entry.composition.num_atoms)
						new_pd_entries.append(entry)
		
		new_pd_entries = self.process_elements_hull(composition, temperature, new_pd_entries)
		
		new_conv_hull = PhaseDiagram(new_pd_entries)
		return new_conv_hull
	
	def process_elements_hull(self, composition, temperature, pd_entries_list):
		for ele in composition:
			if ele in ['Fe', 'Ti', 'Mn', 'Hf', 'Zr'] and self.correction:
				
				for key, value in self.end_member[ele].items():
					name = Composition(ele)
					if key == self.transition_temperatures[ele][1]:
						temp_energy = value - value * temperature / self.transition_temperatures[ele][2]
					else:
						temp_energy = value
					
					pd_entries_list.append(PDEntryLocal(composition=name, energy=temp_energy * name.num_atoms,
														name=f'{name.alphabetical_formula}_{key}'))
			else:
				for key, value in self.end_member[ele].items():
					name = Composition(ele)
					
					pd_entries_list.append(PDEntryLocal(composition=name, energy=value * name.num_atoms,
														name=f'{name.alphabetical_formula}_{key}'))
		
		return pd_entries_list
	
	def find_enthalpy_entropy_composition(self,
										  composition: list[str],
										  lattice: str,
										  mol_ratio: list,
										  temperature: float,
										  ):
		mol_ratio = dict(zip(composition, mol_ratio))
		mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}
		
		if len(mol_ratio.keys()) == 1:
			mix_enthalpy = 0
		else:
			
			mix_enthalpy = self.tm.calc_mutinary_multilattice_mix_Enthalpy(mol_ratio=mol_ratio,
																		   binary_dict=self.data,
																		   end_member_dict=self.end_member,
																		   correction=self.correction,
																		   temperature=temperature,
																		   transition_temperatures=self.transition_temperatures,
																		   model='regular')
			if lattice == 'min':
				mix_enthalpy = min(list(mix_enthalpy.values()))
			else:
				mix_enthalpy = mix_enthalpy[lattice]
		
		entropy = self.tm.calc_configEntropy(mol_ratio)
		return mix_enthalpy, entropy, mol_ratio
	
	@staticmethod
	def _make_PD_entry(mol_ratio: dict,
					   energy: float):
		"""
		Helper function to make a PD entry that can be fed into a Phase Diagram.
		Args:
			mol_ratio: list of molar ratios for the composition
			energy: Energy of the phase in eV/atom

		Returns: A PDEntry object

		"""
		
		name = Composition(Composition(mol_ratio).get_integer_formula_and_factor()[0])
		return PDEntryLocal(composition=name, energy=energy * name.num_atoms)
