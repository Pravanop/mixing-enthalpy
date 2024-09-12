from typing import Union

import numpy as np
import pandas as pd
from scipy.integrate import quad

from calculateEnthalpy.helper_functions.grid_code import create_multinary
from calculateEnthalpy.helper_functions.data_utils import DataUtils


class thermoMaths:
	"""
    Set of functions to calculate thermodynamic properties.
    This is from Zhang et al., A Fast and Robust Method for Predicting the Phase Stability of Refractory Complex
    Concentrated Alloys Using Pairwise Mixing Enthalpy http://dx.doi.org/10.2139/ssrn.4081906 """
	
	def __init__(self, ):
		# self.binary_dict = binary_dict
		self.kb = 8.617333262e-05
	
	@staticmethod
	def calc_pairwiseInteractionParameter(
			mix_enthalpy: float,
			mol_i: float,
			mol_j: float) -> float:
		"""
        Function for the Eq. 1 in Zhang et al. A Fast and Robust Method for Predicting the Phase Stability of Refractory
        Complex Concentrated Alloys Using Pairwise Mixing Enthalpy. http://dx.doi.org/10.2139/ssrn.4081906

                                    omega_ij = delta_H_ij/(x_i*x_j)

        Args:
            mix_enthalpy: binary equimolar mixing enthalpy of elements i, j
            mol_i: mole fraction of element i
            mol_j: mole fraction of element j

        Returns: pairwise interaction parameter in eV/atom

        """
		return np.round(mix_enthalpy / (mol_i * mol_j), 5)
	
	@staticmethod
	def calc_subregular_model_enthalpy(
			mol_fraction: list[float],
			omega1: float,
			omega2: float,
	) -> np.array:
		"""
        Calculate a binary subregular model with cubic fit. Taken from John Cavin's thesis Eq 3.3.

        Args:
            mol_fraction: list of molar fractions for the composition in that order
            omega1: pairwise-interaction parameter 1
            omega2: pairwise-interaction parameter 2

        Returns: mixing enthalpy value for mole fraction in eV/atom

        """
		x_i, x_j = mol_fraction
		return x_i * x_j * (omega1 * x_i + omega2 * x_j)
	
	@staticmethod
	def calc_regular_model_enthalpy(
			mol_fraction: list[float],
			omega: float,
	) -> np.array:
		"""
        Calculate a subregular model with cubic fit. Taken from John Cavin's thesis Eq 3.3.

        Args:
            mol_fraction: list of molar fractions for the composition in that order
            omega: the pairwise-interaction parameter in eV/atom

        Returns: binary regular model enthalpy in eV/atom

        """
		x_i, x_j = mol_fraction
		return x_i * x_j * omega
	
	def calc_configEntropy(self,
						   mol_ratio: dict[str:float]) -> float:
		
		"""
        A simple function to calculate boltzmann configurational entropy.

                        delta_S = -k_b*sum_i=1_n(x_i*ln(x_i))

        Args:
            mol_ratio: dictionary of mole fractions corresponding to each element

        Returns: configurational entropy in eV/atom/K

        """
		assert round(sum(list(mol_ratio.values())),
					 3) == 1  # just in case, should be handled internally, but if non-equimolar is used
		
		return np.round(-self.kb * sum([value * np.log(value) for key, value in mol_ratio.items()]), 7)
	
	@staticmethod
	def calc_gibbs_energy(
			enthalpy: float,
			entropy: float,
			temperature: float
	) -> float:
		"""
        Function to implement G = H - TS

        Args:
            enthalpy: mixing enthalpy in eV/atom
            entropy: configurational entropy in eV/atom
            temperature: desirec temperature in K

        Returns: gibbs free energy in eV/atom

        """
		return np.round(enthalpy - temperature * entropy, 5)
	
	def calc_vib_gibbs_energy(self):
		pass
	
	def calc_debye_function(self,
							x: float = 0) -> float:
		"""

        :param x: upper limit for function integrand
        :return:
        """
		debye = 3 / (x ** 3)
		
		def integrand(t) -> float:
			"""

            :param t:
            :return:
            """
			return (t ** 3) / (np.exp(t) - 1)
		
		integral = quad(integrand, 0, x)[0]
		return debye * integral
	
	def calc_vibrational_energy(self,
								debye_temp: float,
								temp: float
								) -> float:
		"""

        :param debye_temp:
        :param temp:
        :return:
        """
		x = debye_temp / temp
		energy = (9 / 8) * self.kb * debye_temp
		energy = energy + self.kb * temp * (3 * np.log(1 - np.exp(-x)) - self.calc_debye_function(x))
		return energy
	
	def avg_T_melt(self, composition, mol_ratio):
		self.meltT = pd.read_csv(
			"/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/data/input_data/PubChemElements_all.csv").to_numpy()
		self.meltT = dict(zip(self.meltT[:, 0], self.meltT[:, 1]))
		tm = [self.meltT[ind] * mol_ratio[i] for i, ind in enumerate(composition)]
		return int(sum(tm))
	
	def calc_mutinary_multilattice_mix_Enthalpy(self,
												mol_ratio,
												binary_dict,
												end_member_dict,
												model: str = "regular",
												correction: bool = True,
												temperature: float = 0) -> Union[dict, int]:
		
		ele_list = list(mol_ratio.keys())
		if len(ele_list) <= 1:
			return 0
		
		binaries = create_multinary(element_list=ele_list, no_comb=[2])
		mix_enthalpy = {}
		if binaries:
			for idx, binary in binaries.items():
				for idx2, ele_pair in enumerate(binary):
					two_el = ele_pair.split("-")
					mix_enthalpy_values = binary_dict[ele_pair]
					for lattice, enthalpy in mix_enthalpy_values.items():
						if not correction:
							if model == "regular":
								omega_ij = self.calc_pairwiseInteractionParameter(mix_enthalpy=enthalpy,
																				  mol_i=0.5,
																				  mol_j=0.5)
						
						else:
							omega_ij = enthalpy  #biased values are stored as omegas, needs to be standaradized
						
						mol_fraction = [mol_ratio[two_el[0]], mol_ratio[two_el[1]]]
						H_mix = self.calc_regular_model_enthalpy(mol_fraction=mol_fraction,
																 omega=omega_ij)
						
							# H_mix += (temp_energies[two_el[0]] - temp_energies[two_el[1]]) * mol_fraction[0]
							
						if lattice not in mix_enthalpy:
							mix_enthalpy[lattice] = H_mix
						else:
							mix_enthalpy[lattice] += H_mix
		
			
			for lattice, enthalpy in mix_enthalpy.items():
				if correction:
					correction_term = 0
					temp_energies = {}
					for end_member in ele_list:
						if end_member in ['Fe', 'Ti', 'Mn', 'Hf', 'Zr']:
							# print(end_member)
							transition_temperatures = {
								'Fe': ['BCC', 'FCC', 1180],
								'Ti': ['HCP', 'BCC', 1155],
								'Hf': ['HCP', 'BCC', 2016],
								'Zr': ['HCP', 'BCC', 1136],
								'Mn': ['BCC', 'FCC', 1370]
							}
							if lattice == transition_temperatures[end_member][1]:  # if the phase has T transition
								temp_energy = end_member_dict[end_member][lattice] - \
											  end_member_dict[end_member][lattice] * float(temperature) / \
											  transition_temperatures[end_member][2]  # T correction is made
								temp_energies[end_member] = temp_energy
							else:
								temp_energy = end_member_dict[end_member][lattice]
								temp_energies[end_member] = temp_energy
						else:
							temp_energy = end_member_dict[end_member][lattice]
							temp_energies[end_member] = temp_energy
					# biasing correction
					for key, value in temp_energies.items():
						
						correction_term += value * mol_ratio[key]
					
					mix_enthalpy[lattice] += correction_term
				
				else:
					continue
			
			return mix_enthalpy
	
	def find_subset_enthalpies(self,
							   ele_list,
							   binary_dict,
							   model: str = "regular") -> pd.DataFrame:
		# ele_list = list(mol_ratio.keys())
		binaries = create_multinary(element_list=ele_list, no_comb=list(range(2, len(ele_list) + 1)))
		mix_enthalpy = {}
		for idx, dimensionality in binaries.items():
			for idx2, ele in enumerate(dimensionality):
				mol_ratio = [1 / idx] * idx
				comp = ele.split('-')
				composition_dict = dict(zip(comp, mol_ratio))
				composition_dict = {key: val for key, val in composition_dict.items() if val != 0.0}
				mix_enthalpy[ele] = self.calc_mutinary_multilattice_mix_Enthalpy(composition_dict, binary_dict, model)
		df = pd.DataFrame(mix_enthalpy)
		return df.T
	
	def calc_multinary_mixEnthalpy(self,
								   alloy_comp: str,
								   mol_ratio: dict[str:float],
								   **kwargs
								   ) -> float:
		"""
        Generalizing to off equimolar compositions with regular model

        Args:
            alloy_comp: str of the form 'ele1-ele2-ele3'
            mol_ratio: mole fraction for each element in composition in order

        Returns: enthalpy of mixing for multinary compound in ev/atom

        """
		ele_list = alloy_comp.split('-')
		binaries = create_multinary(element_list=ele_list, no_comb=[2])
		mix_enthalpy = 0
		for idx, binary in binaries.items():
			for idx2, ele in enumerate(binary):
				two_el = ele.split("-")
				if self.binary_dict is not None:
					mix_enthalpy_data = DataUtils().extract_binaryEnthalpy(
						binary_dict=self.binary_dict,
						ele_pair=ele
					)
				else:
					if kwargs["binary_dict"] is not None:
						mix_enthalpy_data = DataUtils().extract_binaryEnthalpy(
							binary_dict=kwargs["binary_dict"],
							ele_pair=ele
						)
				omega_ij = self.calc_pairwiseInteractionParameter(mix_enthalpy=mix_enthalpy_data,
																  mol_i=0.5,
																  mol_j=0.5)
				mol_fraction = [mol_ratio[two_el[0]], mol_ratio[two_el[1]]]
				H_mix = self.calc_regular_model_enthalpy(mol_fraction=mol_fraction,
														 omega=omega_ij)
				mix_enthalpy += H_mix
		return np.round(mix_enthalpy, 6)
