"""This is from Zhang et al., A Fast and Robust Method for Predicting the Phase Stability of Refractory Complex
Concentrated Alloys Using Pairwise Mixing Enthalpy http://dx.doi.org/10.2139/ssrn.4081906 """
import numpy as np
from calculateEnthalpy.create_alloy_comp import create_multinary
from calculateEnthalpy.data_utils import extract_binaryEnthalpy
def calc_pairwiseInteractionParameter(mix_enthalpy: float , mol_i: float , mol_j: float) -> float :
	"""
	Function for the Eq. 1 in Zhang et al. A Fast and Robust Method for Predicting the Phase Stability of Refractory
	Complex Concentrated Alloys Using Pairwise Mixing Enthalpy. http://dx.doi.org/10.2139/ssrn.4081906

								omega_ij = delta_H_ij/(x_i*x_j)

	:param mix_enthalpy: binary equimolar mixing enthalpy of elements i, j
	:param mol_i: mole fraction of element i
	:param mol_j: mole fraction of element j
	:return: pairwise interaction parameter
	"""
	return np.round(mix_enthalpy / (mol_i * mol_j) , 7)

def gibbs_energy(enthalpy: float,
                 entropy: float,
                 temperature: float) -> float:
	"""
	Function to implement G = H - TS
	:param enthalpy: the mixing enthalpy
	:param entropy: the configurational entropy
	:param temperature: the desired temperature
	:return: gibbs free energy
	"""
	return np.round(enthalpy - temperature*entropy, 7)

def calc_multinary_mixEnthalpy(
		alloy_comp: str ,
		binary_dict: dict ,
		mol_ratio: dict ,
		) -> float :
	"""
	Function for Eq. 2 in Zhang et al. A Fast and Robust Method for Predicting the Phase Stability of Refractory
	Complex Concentrated Alloys Using Pairwise Mixing Enthalpy. http://dx.doi.org/10.2139/ssrn.4081906

							delta_Hmix = sum_i=1,i!=j_n(omega_ij*x_i*x_j)

	:param binary_dict:
	:param alloy_comp: str of the form 'ele1-ele2-ele3'
	:param mol_ratio: mole fraction for each element in composition in order!
	:return: enthalpy of mixing for multinary compound in ev/atom
	"""
	ele_list = alloy_comp.split('-')
	binaries = create_multinary(element_list = ele_list , no_comb = [2])[0]
	mix_enthalpy = 0
	for idx , binary in enumerate(binaries) :
		binary_mole_ratio = 0.5
		delH_ij = extract_binaryEnthalpy(
				binary_dict = binary_dict ,
				ele_pair = binary
				)
		omega_ij = calc_pairwiseInteractionParameter(
				mix_enthalpy = delH_ij ,
				mol_i = binary_mole_ratio ,
				mol_j = 1 - binary_mole_ratio
				)
		two_el = binary.split("-")
		mix_enthalpy += omega_ij * mol_ratio[two_el[0]] * mol_ratio[two_el[1]]  # the implementation of the main formula
	
	return np.round(mix_enthalpy , 7)

def calc_configEntropy(mol_ratio: dict) -> float :
	"""
	A simple function to calculate boltzmann configurational entropy.

					delta_S = -k_b*sum_i=1_n(x_i*ln(x_i))

	:param mol_ratio: list of floats of mole fraction
	:return: boltzmann configurational entropy
	"""
	assert sum(list(mol_ratio.values())) == 1  # just in case, should be handled internally, but if non-equimolar is used
	k_b = 8.617333262e-05
	return np.round(-k_b * sum([value * np.log(value) for key , value in mol_ratio.items()]) , 7)  # implementation of
# main formula