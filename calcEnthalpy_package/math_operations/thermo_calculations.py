from typing import Union

import numpy as np
import pandas as pd
from scipy.integrate import quad
from calcEnthalpy_package.grids_and_combinations.combination_generation import MultinaryCombinations  # Updated import


class ThermoMaths:
	"""
	A set of functions to calculate thermodynamic properties such as enthalpy, entropy, Gibbs free energy, and vibrational energy.

	Based on Zhang et al., "A Fast and Robust Method for Predicting the Phase Stability of Refractory Complex
	Concentrated Alloys Using Pairwise Mixing Enthalpy". http://dx.doi.org/10.2139/ssrn.4081906
	"""
	
	def __init__(self):
		"""
		Initializes constants like the Boltzmann constant in eV/K.
		"""
		self.kb = 8.617e-05  # Boltzmann constant in eV/K
		meltT_path = pd.read_csv(
			"../database/PubChemElements_all.csv").to_numpy()
		self.meltT = dict(zip(meltT_path[:, 0], meltT_path[:, 1]))
		
	@staticmethod
	def calc_pairwiseInteractionParameter(
			mix_enthalpy: float,
			mol_i: float,
			mol_j: float) -> float:
		"""
		Calculates the pairwise interaction parameter (omega_ij) from the binary mixing enthalpy.

		Parameters
		----------
		mix_enthalpy : float
			Binary equimolar mixing enthalpy of elements i and j.
		mol_i : float
			Mole fraction of element i.
		mol_j : float
			Mole fraction of element j.

		Returns
		-------
		float
			Pairwise interaction parameter in eV/atom.
		"""
		return mix_enthalpy / (mol_i * mol_j)
	
	def avg_T_melt(self, composition, mol_ratio):
		if isinstance(composition, str):
			tm = self.meltT[composition]
			return tm
		else:
			tm = [self.meltT[ind] * mol_ratio[i] for i, ind in enumerate(composition)]
			return int(sum(tm))
	
	@staticmethod
	def calc_subregular_model_enthalpy(
			mol_fraction: list[float],
			omega1: float,
			omega2: float) -> np.array:
		"""
		Calculates mixing enthalpy using a binary subregular model with cubic fit.

		Parameters
		----------
		mol_fraction : list of float
			List of molar fractions for the composition.
		omega1 : float
			Pairwise-interaction parameter 1.
		omega2 : float
			Pairwise-interaction parameter 2.

		Returns
		-------
		np.array
			Mixing enthalpy value in eV/atom.
		"""
		x_i, x_j = mol_fraction
		return x_i * x_j * (omega1 * x_i + omega2 * x_j)
	
	@staticmethod
	def calc_regular_model_enthalpy(
			mol_fraction: list[float],
			omega: float) -> np.array:
		"""
		Calculates mixing enthalpy using a regular model.

		Parameters
		----------
		mol_fraction : list of float
			List of molar fractions for the composition.
		omega : float
			Pairwise-interaction parameter in eV/atom.

		Returns
		-------
		np.array
			Binary regular model enthalpy in eV/atom.
		"""
		x_i, x_j = mol_fraction
		return x_i * x_j * omega
	
	def calc_configEntropy(self, mol_ratio: dict[str, float]) -> float:
		"""
		Calculates Boltzmann configurational entropy.

		Parameters
		----------
		mol_ratio : dict
			Dictionary of mole fractions corresponding to each element.

		Returns
		-------
		float
			Configurational entropy in eV/atom/K.
		"""
		assert round(sum(mol_ratio.values()), 3) == 1, "Mole fractions must sum to 1."
		
		entropy = -self.kb * sum([value * np.log(value) for value in mol_ratio.values()])
		return entropy
	
	@staticmethod
	def calc_gibbs_energy(enthalpy: float, entropy: float, temperature: float) -> float:
		"""
		Calculates Gibbs free energy using G = H - TS.

		Parameters
		----------
		enthalpy : float
			Mixing enthalpy in eV/atom.
		entropy : float
			Configurational entropy in eV/atom/K.
		temperature : float
			Temperature in K.

		Returns
		-------
		float
			Gibbs free energy in eV/atom.
		"""
		return enthalpy - temperature * entropy
	
	def calc_debye_function(self, x: float = 0) -> float:
		"""
		Calculates the Debye function for vibrational energy calculation.

		Parameters
		----------
		x : float
			Upper limit for function integrand.

		Returns
		-------
		float
			Value of the Debye function.
		"""
		debye = 3 / (x ** 3)
		
		def integrand(t) -> float:
			"""Integrand function for Debye integral."""
			return (t ** 3) / (np.exp(t) - 1)
		
		integral = quad(integrand, 0, x)[0]
		return debye * integral
	
	def calc_vibrational_energy(self, debye_temp: float, temp: float) -> float:
		"""
		Calculates vibrational energy using the Debye model.

		Parameters
		----------
		debye_temp : float
			Debye temperature in K.
		temp : float
			Temperature in K.

		Returns
		-------
		float
			Vibrational energy in eV.
		"""
		x = debye_temp / temp
		energy = (9 / 8) * self.kb * debye_temp
		energy += self.kb * temp * (3 * np.log(1 - np.exp(-x)) - self.calc_debye_function(x))
		return energy
	
	def calc_mutinary_multilattice_mix_Enthalpy(self,
												mol_ratio: dict[str:float],
												binary_dict: dict[str:float],
												end_member_dict: dict[str:float],
												transition_temperatures: dict[str, list],
												correction: bool,
												temperature: float,
												model: str) -> Union[dict, int]:
		
		ele_list = list(mol_ratio.keys())
		if len(ele_list) <= 1:
			return 0
		
		binaries = MultinaryCombinations.create_multinary(element_list=ele_list, no_comb=[2])
		mix_enthalpy = {}
		
		if binaries:
			
			for binary in binaries.values():
				for ele_pair in binary:
					two_el = ele_pair.split("-")
					mix_enthalpy_values = binary_dict[ele_pair]
					
					# Precompute mol_fraction once for the pair
					mol_fraction = [mol_ratio[two_el[0]], mol_ratio[two_el[1]]]
					
					for lattice, enthalpy in mix_enthalpy_values.items():
						# Calculate interaction parameter or use pre-biased values
						if not correction:
							if model == "regular":
								omega_ij = self.calc_pairwiseInteractionParameter(
									mix_enthalpy=enthalpy, mol_i=0.5, mol_j=0.5
								)
							else:
								omega_ij = self.calc_pairwiseInteractionParameter(
									mix_enthalpy=enthalpy, mol_i=0.5, mol_j=0.5
								)
						else:
							omega_ij = enthalpy  # biased values
						
						# Calculate enthalpy using the regular model
						H_mix = self.calc_regular_model_enthalpy(
							mol_fraction=mol_fraction, omega=omega_ij
						)
						
						if correction:
							# Biasing correction for transition metals
							temp_energies = {}
							for end_member in two_el:
								if end_member in transition_temperatures:
									# Apply temperature correction for specific metals
									transition_data = transition_temperatures[end_member]
									if lattice == transition_data[1]:  # Phase with T transition
										temp_energy = end_member_dict[end_member][lattice] - \
													  (end_member_dict[end_member][lattice] *
													   float(temperature) / transition_data[2])
									else:
										temp_energy = end_member_dict[end_member][lattice]
								else:
									temp_energy = end_member_dict[end_member][lattice]
								
								# Store temp energy for each element
								temp_energies[end_member] = temp_energy
							
							# Apply the biasing correction to H_mix
							H_mix += temp_energies[two_el[1]] * mol_fraction[1] + \
									 temp_energies[two_el[0]] * mol_fraction[0]
						
						# Add enthalpy contribution to the mix_enthalpy dictionary
						if lattice not in mix_enthalpy:
							mix_enthalpy[lattice] = H_mix
						else:
							mix_enthalpy[lattice] += H_mix
			
			return mix_enthalpy
