from typing import Union, Tuple, Any
import numpy as np
import pandas as pd
from numpy import ndarray, dtype, signedinteger
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Element
from tqdm import tqdm

from calcEnthalpy_package.grids_and_combinations.combination_generation import MultinaryCombinations
from calcEnthalpy_package.grids_and_combinations.grid_creation import CompositionGrid
from calcEnthalpy_package.math_operations.thermo_calculations import ThermoMaths
from calcEnthalpy_package.phase_diagram.convex_hull import ConvexHull

from calcEnthalpy_package.data_processing.intermetallic_extraction import IntermetallicExtractions


class GridIterator:
	
	def __init__(self,
				 grid_size: int,
				 tm: ThermoMaths,
				 data: dict,
				 end_member: dict,
				 flags: dict,
				 api_key: str
				 ):
		
		self.grid_size = grid_size
		self.end_member = end_member
		self.im_flag = flags['im_flag']
		self.equi_flag = flags['equi_flag']
		self.correction = flags['correction']
		self.tm = tm
		self.api_key = api_key
		self.data = data
		self.convex_hull = ConvexHull(
			end_member=end_member,
			grid_size=grid_size,
			flags=flags,
			tm=tm,
			api_key=api_key,
			data=data
		)
		
		self.temp_grid = None
	
	def upper_limit(self, composition, mol_ratio, phase_flag) -> None:
		avg_tm = self.tm.avg_T_melt(composition, mol_ratio)
		if phase_flag:
			self.temp_grid = np.arange(0, avg_tm + 200, 20, dtype=np.float64)
		else:
			self.temp_grid = np.arange(0, avg_tm + 100, 100, dtype=np.float64)
	
	def temp_iterator(self,
					  composition: list,
					  temp_grid: list) -> dict[float:PhaseDiagram]:
		
		n_alloy = len(composition)
		all_combs = MultinaryCombinations.create_multinary(element_list=composition,
														   no_comb=list(range(2, n_alloy + 1)))
		
		im_list = []
		for dimensionality, alloy_list in all_combs.items():
			if self.im_flag:
				im_list += IntermetallicExtractions.get_MP_intermetallic(alloy_list, api_key=self.api_key)
		
		PD_temp_comp_dict = {}
		for idx, temp in enumerate(tqdm(temp_grid, desc="Running Temperature")):
			PD_temp_comp_dict[temp] = self.convex_hull.make_convex_hull(temperature=temp,
																		composition=composition,
																		batch_tag=True,
																		im=im_list)
		
		return PD_temp_comp_dict
	
	def misc_temperature_across_grid(self,
									 composition: list[str],
									 mol_grid_size: int,
									 lattice: str,
									 phase_flag: bool = False,
									 **kwargs
									 ) -> tuple[ndarray, ndarray]:
		
		n = len(composition)
		
		if isinstance(mol_grid_size, int):
			mol_grid = CompositionGrid.create_mol_grid(n, mol_grid_size)
		else:
			mol_grid = mol_grid_size
		
		temp_ele = []
		for i in composition:
			temp_ele.append(self.tm.avg_T_melt(i, mol_grid))
		
		
		
		if 'conv_hull' in kwargs:
			conv_hull_dict = kwargs['conv_hull']
			temp_grid = kwargs['temp_grid']
		else:
			if phase_flag:
				temp_grid = list(np.arange(0, max(temp_ele) + 100, 20, dtype=np.float64))
			else:
				temp_grid = list(np.arange(0, max(temp_ele) + 100, 100, dtype=np.float64))
			conv_hull_dict = self.temp_iterator(composition=composition,
												temp_grid=temp_grid)
		
		stables = np.ones_like(mol_grid[:, 0]) * -1
		
		done_idx = []
		for idx_temp, temp in enumerate(temp_grid):
			
			for idx, mol in enumerate(mol_grid):
				if idx in done_idx:
					continue
				
				mix_enthalpy, entropy, mol_ratio = self.convex_hull.find_enthalpy_entropy_composition(
					mol_ratio=list(mol),
					composition=composition,
					temperature=temp,
					lattice=lattice,
				)
				is_stable = self.get_decomp_and_e_hull(mol_ratio=mol_ratio,
													   temp=temp,
													   conv_hull=conv_hull_dict[temp],
													   mix_enthalpy=mix_enthalpy,
													   entropy=entropy)
				
				if is_stable is not None:
					# Using a tolerance of 1 meV/atom for computing energy above the hull
					if np.isclose(is_stable[1], 0.0, atol=1e-3):
						# print(mol, is_stable)
						stables[idx] = temp
						done_idx.append(idx)
					else:
						continue
				else:
					stables[idx] = temp
					done_idx.append(idx)
		
		return mol_grid, stables
	
	def e_hull_across_grid(self,
						   composition: list[str],
						   mol_grid_size: int,
						   lattice: str,
						   single_temp_flag: bool = False,
						   temp_gradation: int = 300,
						   ) -> tuple[ndarray, ndarray, ndarray]:
		
		n = len(composition)
		if isinstance(mol_grid_size, int):
			mol_grid = CompositionGrid.create_mol_grid(n, mol_grid_size)
		else:
			mol_grid = mol_grid_size
		
		temp_ele = []
		for i in composition:
			temp_ele.append(self.tm.avg_T_melt(i, mol_grid))
		
		if single_temp_flag:
			temp_grid = [temp_gradation]
		else:
			temp_grid = list(np.arange(200, max(temp_ele) + 100, temp_gradation, dtype=np.float64))
		
		conv_hull_dict = self.temp_iterator(composition=composition,
											temp_grid=temp_grid)
		
		stables = np.ones((len(mol_grid), len(temp_grid)))
		
		for idx_temp, temp in enumerate(temp_grid):
			
			for idx, mol in enumerate(mol_grid):
				
				mix_enthalpy, entropy, mol_ratio = self.convex_hull.find_enthalpy_entropy_composition(
					mol_ratio=list(mol),
					composition=composition,
					temperature=temp,
					lattice=lattice,
				)
				is_stable = self.get_decomp_and_e_hull(mol_ratio=mol_ratio,
													   temp=temp,
													   conv_hull=conv_hull_dict[temp],
													   mix_enthalpy=mix_enthalpy,
													   entropy=entropy)
				
				if is_stable is not None:
					# Using a tolerance of 1 meV/atom for computing energy above the hull
					if np.isclose(is_stable[1], 0.0, atol=1e-3):
						stables[idx][idx_temp] = 0.0
					else:
						stables[idx][idx_temp] = is_stable[1]
				else:
					stables[idx][idx_temp] = 0.0
		
		return mol_grid, stables, temp_grid
	
	def uni_molar_misc_temperature(self,
								   mol_ratio: list,
								   composition: list[str],
								   lattice: str,
								   phase_flag: bool = False,
								   batch_tag: bool = False,
								   **kwargs) -> Union[float, tuple[float, int]]:
		
		self.upper_limit(composition, mol_ratio, phase_flag)
		
		if not batch_tag:
			n_alloy = len(composition)
			all_combs = MultinaryCombinations.create_multinary(element_list=composition,
															   no_comb=list(range(2, n_alloy + 1)))
			
			im_list = []
			for dimensionality, alloy_list in all_combs.items():
				if self.im_flag:
					im_list += IntermetallicExtractions.get_MP_intermetallic(alloy_list,
																			 api_key=self.api_key)
		else:
			im_list = []
			if self.im_flag:
				im_list = kwargs.get('im', [])
		
		if "conv_hull" not in kwargs:
			conv_hull0 = self.convex_hull.make_convex_hull(
				temperature=0,
				composition=composition,
				batch_tag=True,
				im=im_list)
		else:
			conv_hull0 = kwargs['conv_hull']
		
		for idx, temperature in enumerate(self.temp_grid):
			mix_enthalpy, entropy, mol_ratio = self.convex_hull.find_enthalpy_entropy_composition(
				mol_ratio=mol_ratio,
				composition=composition,
				temperature=float(temperature),
				lattice=lattice,
			)
			
			conv_hull = self.convex_hull.update_convex_hull(
				composition=composition,
				conv_hull=conv_hull0,
				temperature=temperature)
			
			is_stable = self.get_decomp_and_e_hull(mol_ratio=mol_ratio,
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
		
		return float(self.temp_grid[-1]), 0  # to demarcate melting point
	
	def find_decomp_products(self,
							 mol_ratio: list,
							 composition: list[str],
							 temperature: float,
							 lattice: str,
							 batch_tag: bool = False,
							 **kwargs) -> Union[tuple, int]:
		
		full_comp = composition
		
		mix_enthalpy, entropy, mol_ratio = self.convex_hull.find_enthalpy_entropy_composition(
			mol_ratio=mol_ratio,
			composition=composition,
			temperature=temperature,
			lattice=lattice,
		)
		
		if batch_tag and self.im_flag:
			im_list = kwargs.get('im', [])
		else:
			raise 'Provide Intermetallics'
		
		entropy = self.tm.calc_configEntropy(mol_ratio)
		composition = list(mol_ratio.keys())
		exclude = list(set(full_comp).difference(set(composition)))
		im_list_exclude = []
		if exclude:
			for j in exclude:
				for idx, i in enumerate(im_list):
					if Element(j) in i.elements:
						im_list_exclude.append(i)
		
		im_list = list(set(im_list).difference(im_list_exclude))
		if batch_tag:
			conv_hull = self.convex_hull.make_convex_hull(temperature=float(temperature),
														  composition=composition,
														  batch_tag=True,
														  im=im_list
														  )
		else:
			conv_hull = self.convex_hull.make_convex_hull(temperature=float(temperature),
														  composition=composition,
														  )
		
		is_stable = self.get_decomp_and_e_hull(mol_ratio=mol_ratio,
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
	
	def get_decomp_and_e_hull(self,
							  mol_ratio: dict,
							  temp: float,
							  conv_hull: PhaseDiagram,
							  mix_enthalpy: float,
							  entropy: float) -> Union[tuple[dict, float], None]:
		
		gibbs = self.tm.calc_gibbs_energy(enthalpy=mix_enthalpy,
										  entropy=entropy,
										  temperature=temp)
		pdEntry = self.convex_hull._make_PD_entry(energy=gibbs,
												  mol_ratio=mol_ratio)
		try:
			
			answer = conv_hull.get_decomp_and_e_above_hull(entry=pdEntry,
														   check_stable=False)
			return answer
		
		except ValueError as e:
			return None
