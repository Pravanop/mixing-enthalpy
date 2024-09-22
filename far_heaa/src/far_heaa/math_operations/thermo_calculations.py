from typing import Union, List, Dict
import numpy as np
import pandas as pd
from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations


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
        return -self.kb * sum([value * np.log(value) for value in mol_ratio.values()])

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

        binaries = MultinaryCombinations.create_multinary(
            element_list=ele_list, no_comb=[2]
        )
        mix_enthalpy = {}

        if binaries:
            for binary in binaries.values():
                for ele_pair in binary:
                    two_el = ele_pair.split("-")
                    mix_enthalpy_values = binary_dict[ele_pair]

                    mol_fraction = [mol_ratio[two_el[0]], mol_ratio[two_el[1]]]

                    for lattice, enthalpy in mix_enthalpy_values.items():
                        if not correction:
                            omega_ij = self.calc_pairwise_interaction_parameter(
                                mix_enthalpy=enthalpy, mol_i=0.5, mol_j=0.5
                            )
                        else:
                            omega_ij = enthalpy  # biased values

                        H_mix = self.calc_regular_model_enthalpy(
                            mol_fraction=mol_fraction, omega=omega_ij
                        )

                        if correction:
                            temp_energies = {}
                            for end_member in two_el:
                                if end_member in transition_temperatures:
                                    transition_data = transition_temperatures[
                                        end_member
                                    ]
                                    if lattice == transition_data[1]:
                                        temp_energy = end_member_dict[end_member][
                                            lattice
                                        ] - (
                                            end_member_dict[end_member][lattice]
                                            * float(temperature)
                                            / transition_data[2]
                                        )
                                    else:
                                        temp_energy = end_member_dict[end_member][
                                            lattice
                                        ]
                                else:
                                    temp_energy = end_member_dict[end_member][lattice]

                                temp_energies[end_member] = temp_energy

                            H_mix += (
                                temp_energies[two_el[1]] * mol_fraction[1]
                                + temp_energies[two_el[0]] * mol_fraction[0]
                            )

                        if lattice not in mix_enthalpy:
                            mix_enthalpy[lattice] = H_mix
                        else:
                            mix_enthalpy[lattice] += H_mix

            return mix_enthalpy
