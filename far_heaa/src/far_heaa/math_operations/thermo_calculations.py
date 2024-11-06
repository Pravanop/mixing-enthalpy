from typing import Union, List, Dict
import numpy as np
import pandas as pd
from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from sympy import symbols
from sympy.functions.elementary.exponential import log
from sympy.parsing.sympy_parser import parse_expr
from sympy.matrices.dense import Matrix


class ThermoMaths:
    """
    A set of functions to calculate thermodynamic properties such as enthalpy, entropy, Gibbs free energy,
    and vibrational energy.

    Based on Zhang et al., "A Fast and Robust Method for Predicting the Phase Stability of Refractory Complex
    Concentrated Alloys Using Pairwise Mixing Enthalpy". http://dx.doi.org/10.2139/ssrn.4081906
    """
    
    def __init__(self):
        """
        Initializes the ThermoMaths class, setting the Boltzmann constant and loading the melting temperature
        database.
        """
        self.kb = 8.617e-05  # Boltzmann constant in eV/K
        melt_T_path = pd.read_csv("../database/PubChemElements_all.csv").to_numpy()
        self.melt_T = dict(zip(melt_T_path[:, 0], melt_T_path[:, 1]))
    
    @staticmethod
    def calc_pairwise_interaction_parameter(
            mix_enthalpy: float, mol_i: float, mol_j: float
    ) -> float:
        """
        Calculates the pairwise interaction parameter (omega_ij) from the binary mixing enthalpy.

        Args:
                mix_enthalpy (float): Binary equimolar mixing enthalpy of elements i and j.
                mol_i (float): Mole fraction of element i.
                mol_j (float): Mole fraction of element j.

        Returns:
                float: Pairwise interaction parameter in eV/atom.
        """
        return mix_enthalpy / (mol_i * mol_j)
    
    def create_xsyms(self, n):
        '''returns a string of symbols: (x1,x2,...,xn-1)
        Contributor: John Cavin'''
        s = ''
        for i in range(n - 1):
            s = s + 'x{} '.format(int(i + 1))
        return symbols(s)
    
    def x_N(self, syms):
        '''input: n-1 syms
        returns 1-x1-x2-...x_(n-1)
        Contributor: John Cavin'''
        xn = parse_expr('1')
        # print(syms)
        for s in syms:
            xn = xn - s
        return xn
    
    def create_Ssym(self, xsyms):
        out = parse_expr('0')
        for i in xsyms:
            out += - tm.kb * i * log(i + 1e-4)
        
        return out
    
    def create_Hsym(self, xsyms, composition, data, lattice):
        binaries = list(MultinaryCombinations.create_multinary(composition, no_comb=[2]).values())[0]
        out = parse_expr('0')
        for binary in binaries:
            omega = data[binary][lattice]
            two_eles = binary.split('-')
            out += omega * xsyms[composition.index(two_eles[0])] * xsyms[composition.index(two_eles[1])]
        
        return out
    
    def hessian(f, syms):
        '''compute the hessian of a matrix of symbols
        Contributor: John Cavin'''
        out = Matrix([[f.diff(x).diff(y) for x in syms] for y in syms])
        return out
    
    def G_sym(self, H_sym, S_sym, T):
        return H_sym - T * S_sym
    
    def find_eigenvalue(self, H_num):
        H_num = np.array(H_num).astype(np.float64)
        w, _ = np.linalg.eig(H_num)
        return w
    
    def avg_T_melt(
            self, composition: Union[str, List[str]], mol_ratio: List[float]
    ) -> float:
        """
        Calculates the average melting temperature for a given alloy composition.

        Args:
                composition (Union[str, List[str]]): The element composition of the alloy.
                mol_ratio (List[float]): The mole fractions of the elements.

        Returns:
                float: The average melting temperature.
        """
        if isinstance(composition, str):
            return self.melt_T[composition]
        
        tm = [self.melt_T[ind] * mol_ratio[i] for i, ind in enumerate(composition)]
        return int(sum(tm))
    
    @staticmethod
    def calc_subregular_model_enthalpy(
            mol_fraction: List[float], omega1: float, omega2: float
    ) -> np.array:
        """
        Calculates the enthalpy using the subregular model.

        Args:
                mol_fraction (List[float]): Mole fractions of the components.
                omega1 (float): The first interaction parameter.
                omega2 (float): The second interaction parameter.

        Returns:
                np.array: The calculated enthalpy.
        """
        x_i, x_j = mol_fraction
        return x_i * x_j * (omega1 * x_i + omega2 * x_j)
    
    @staticmethod
    def calc_regular_model_enthalpy(
            mol_fraction: List[float], omega: float
    ) -> np.array:
        """
        Calculates the enthalpy using the regular model.

        Args:
                mol_fraction (List[float]): Mole fractions of the components.
                omega (float): The interaction parameter.

        Returns:
                np.array: The calculated enthalpy.
        """
        x_i, x_j = mol_fraction
        return x_i * x_j * omega
    
    def calc_config_entropy(self, mol_ratio: Dict[str, float]) -> float:
        """
        Calculates the configurational entropy for a given mole ratio.

        Args:
                mol_ratio (Dict[str, float]): A dictionary of mole fractions for each element.

        Returns:
                float: The configurational entropy in eV/K.

        Raises:
                AssertionError: If the mole fractions do not sum to 1.
        """
        assert round(sum(mol_ratio.values()), 3) == 1, "Mole fractions must sum to 1."
        
        return -self.kb * sum([value * np.log(value + 1e-7) for value in mol_ratio.values()])
    
    @staticmethod
    def calc_gibbs_energy(enthalpy: float, entropy: float, temperature: float) -> float:
        """
        Calculates the Gibbs free energy.

        Args:
                enthalpy (float): The enthalpy value.
                entropy (float): The entropy value.
                temperature (float): The temperature in Kelvin.

        Returns:
                float: The Gibbs free energy.
        """
        return enthalpy - temperature * entropy
    
    def calc_mutinary_multilattice_mix_enthalpy(
            self,
            mol_ratio: Dict[str, float],
            binary_dict: Dict[str, Dict[str, float]],
            end_member_dict: Dict[str, Dict[str, float]],
            transition_temperatures: Dict[str, List[Union[str, float]]],
            correction: bool,
            temperature: float,
            model: str,
    ) -> Union[Dict[str, float], int]:
        """
        Calculates the multinary multilattice mixing enthalpy for an alloy system.

        Args:
                mol_ratio (Dict[str, float]): Mole ratios of the components.
                binary_dict (Dict[str, Dict[str, float]]): Dictionary of binary mixing enthalpies.
                end_member_dict (Dict[str, Dict[str, float]]): Dictionary of end member energies for each phase.
                transition_temperatures (Dict[str, List[Union[str, float]]]): Transition temperatures and phases.
                correction (bool): Whether to apply a correction for transition metals.
                temperature (float): The temperature in Kelvin.
                model (str): The thermodynamic model to use ("regular" or "subregular").

        Returns:
                Union[Dict[str, float], int]: A dictionary of mixing enthalpies for each lattice or 0 if
                                                                          the number of elements is less than 2.
        """
        ele_list = list(mol_ratio.keys())
        if len(ele_list) <= 1:
            return 0
        
        binaries = list(MultinaryCombinations.create_multinary(
            element_list=ele_list, no_comb=[2]
        ).values())[0]
        mix_enthalpy = {}
        
        for ele_pair in binaries:
            two_el = ele_pair.split("-")
            mix_enthalpy_values = binary_dict[ele_pair]
            
            mol_fraction = [mol_ratio[two_el[0]], mol_ratio[two_el[1]]]
            
            for lattice, enthalpy in mix_enthalpy_values.items():
                omega_ij = enthalpy if correction else self.calc_pairwise_interaction_parameter(enthalpy, 0.5,
                                                                                                0.5)
                
                H_mix = (self.calc_regular_model_enthalpy(mol_fraction, omega_ij)
                         if isinstance(omega_ij, float)
                         else self.calc_subregular_model_enthalpy(mol_fraction, omega_ij[1], omega_ij[0]))
                
                if correction:
                    temp_energies = self.temp_correction(end_member_dict,
                                                         two_el,
                                                         lattice,
                                                         transition_temperatures,
                                                         temperature)
                    
                    H_mix += sum(temp_energies[el] * mol_fraction[i] for i, el in enumerate(two_el))
                
                mix_enthalpy[lattice] = mix_enthalpy.get(lattice, 0) + H_mix
            
            return mix_enthalpy
    
    def temp_correction(self,
                        end_member_dict: Dict[str, Dict[str, float]],
                        two_el: List[str],
                        lattice: str,
                        transition_temperatures: Dict[str, List[Union[str, float]]],
                        temperature: float) -> Dict[str, float]:
        temp_energies = {}
        for end_member in two_el:
            temp_energy = end_member_dict[end_member][lattice]
            if end_member in transition_temperatures:
                transition_data = transition_temperatures[end_member]
                if lattice == transition_data[1]:
                    temp_energy -= temp_energy * float(temperature) / transition_data[2]
            
            temp_energies[end_member] = temp_energy
            return temp_energies
    
    def calc_multinary_multilattice_mix_enthalpy_sym(self,
                                                     composition: List[str],
                                                     data: Dict[str, Dict[str, float]],
                                                     correction: bool) -> Union[Dict[str, float], int]:
        
        xsym = list(self.create_xsyms(len(composition)))
        xsyms = xsym + [self.x_N(self.create_xsyms(len(composition)))]
        mix_enthalpy_sym = {}
        for binary in binaries:
            omega_pair = data[binary]
            for lattice in omega_pair:
                omega = omega_pair[lattice]
                two_eles = binary.split('-')
                out = omega * xsyms[composition.index(two_eles[0])] * xsyms[composition.index(two_eles[1])]
                if correction:
                    temp_energies = self.temp_correction(end_member_dict,
                                                         two_el,
                                                         lattice,
                                                         transition_temperatures,
                                                         temperature)
                    
                    out += parse_expr(str(sum(temp_energies[el] * mol_fraction[i] for i, el in enumerate(two_el))))
                
                mix_enthalpy_sym[lattice] = mix_enthalpy.get(lattice, parse_expr('0')) + out
        
        return mix_enthalpy_sym
