from typing import List, Dict, Union, Any

import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Composition

from far_heaa.data_processing.intermetallic_extraction import IntermetallicExtractions
from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.grids_and_combinations.grid_creation import CompositionGrid
from far_heaa.math_operations.thermo_calculations import ThermoMaths
from far_heaa.phase_diagram.pdEntry_local import PDEntryLocal


class ConvexHull:
    """
    A class to construct the convex hull of an alloy system, calculate thermodynamic properties,
    and update the hull with increasing temperature.

    Args:
            end_member (dict): Dictionary containing the end member energies for each element and lattice.
            data (dict): Dictionary containing binary mixing enthalpies.
            flags (dict): Dictionary of flags that determine behavior (e.g., correction, im_flag, equi_flag).
            tm (ThermoMaths): An instance of the ThermoMaths class for thermodynamic calculations.
            api_key (str): API key for accessing external data sources.
            grid_size (int): The size of the molar grid for grid creation.
    """

    def __init__(
        self,
        end_member: Dict,
        data: Dict,
        flags: Dict,
        tm: ThermoMaths,
        api_key: str,
        grid_size: int,
    ):
        """
        Initialize the ConvexHull class with the necessary parameters and data.
        """
        self.transition_temperatures = {
            "Fe": ["BCC", "FCC", 1180],
            "Ti": ["HCP", "BCC", 1155],
            "Hf": ["HCP", "BCC", 2016],
            "Zr": ["HCP", "BCC", 1136],
            "Mn": ["BCC", "FCC", 1370],
        }
        self.end_member = end_member
        self.correction = flags["correction"]
        self.data = data
        self.tm = tm
        self.im_flag = flags["im_flag"]
        self.equi_flag = flags["equi_flag"]
        self.api_key = api_key
        self.grid_size = grid_size

    def make_convex_hull(
        self,
        composition: List[str],
        temperature: float,
        batch_tag: bool = False,
        **kwargs,
    ) -> PhaseDiagram:
        """
        Constructs the convex hull of an alloy system based on the provided composition and temperature.

        Args:
                composition (List[str]): A list of elements in the alloy.
                temperature (float): The temperature in Kelvin for the calculations.
                batch_tag (bool, optional): If True, batch processing of intermetallic phases is enabled. Default is False.
                **kwargs: Additional arguments, such as pre-extracted intermetallic data.

        Returns:
                PhaseDiagram: A PhaseDiagram object representing the convex hull.

        Example::
        
                convex_hull = ConvexHull.make_convex_hull(composition=['Fe', 'Ni'], temperature=1200.0)
        """
        pd_entries_list = []

        if batch_tag and self.im_flag and "im" in kwargs:
            pd_entries_list += kwargs["im"]

        n_alloy = len(composition)
        all_combs = MultinaryCombinations.create_multinary(
            element_list=composition, no_comb=list(range(2, n_alloy + 1))
        )

        for dimensionality, alloy_list in all_combs.items():
            if not batch_tag and self.im_flag:
                pd_entries_list += IntermetallicExtractions.get_MP_intermetallic(
                    alloy_list=alloy_list, api_key=self.api_key
                )
            if self.equi_flag:
                mol_grid = [[1 / dimensionality] * dimensionality]
            else:
                mol_grid = CompositionGrid.create_mol_grid(
                    int(dimensionality), self.grid_size
                )

            for alloy in alloy_list:
                alloy_list_split = alloy.split("-")

                for mol_frac in mol_grid:
                    mol_ratio = dict(zip(alloy_list_split, mol_frac))
                    mol_ratio = {
                        key: val for key, val in mol_ratio.items() if val != 0.0
                    }

                    if len(mol_ratio.keys()) == 1:
                        continue

                    mix_enthalpy = self.tm.calc_mutinary_multilattice_mix_enthalpy(
                        mol_ratio=mol_ratio,
                        binary_dict=self.data,
                        end_member_dict=self.end_member,
                        correction=self.correction,
                        temperature=temperature,
                        transition_temperatures=self.transition_temperatures,
                        model="regular",
                    )
                    config_entropy = self.tm.calc_config_entropy(mol_ratio)

                    if isinstance(mix_enthalpy, dict):
                        for key, enthalpy in mix_enthalpy.items():
                            name = Composition(
                                Composition(mol_ratio).get_integer_formula_and_factor()[
                                    0
                                ]
                            )
                            pd_entries_list.append(
                                PDEntryLocal(
                                    composition=name,
                                    energy=self.tm.calc_gibbs_energy(
                                        enthalpy, config_entropy, temperature
                                    )
                                    * name.num_atoms,
                                    name=f"{name.alphabetical_formula}_{key}",
                                )
                            )

        pd_entries_list = self.process_elements_hull(
            composition, temperature, pd_entries_list
        )
        conv_hull = PhaseDiagram(pd_entries_list)
        return conv_hull

    def update_convex_hull(
        self, composition: List[str], conv_hull: PhaseDiagram, temperature: float
    ) -> PhaseDiagram:
        """
        Updates the convex hull by recalculating the energies of the entries based on the current temperature.

        Args:
                composition (List[str]): A list of elements in the alloy system.
                conv_hull (PhaseDiagram): The current PhaseDiagram object representing the convex hull.
                temperature (float): The temperature in K for the calculations.

        Returns:
                PhaseDiagram: The updated PhaseDiagram object.
        """
        pd_entries = conv_hull.all_entries
        new_pd_entries = []

        for entry in pd_entries:
            lattice = entry.name.split("_")[-1]
            if "MP" in entry.name:
                new_pd_entries.append(entry)
            else:
                mol_ratio = {
                    el.symbol: np.round(entry.composition.get_atomic_fraction(el), 4)
                    for el in entry.composition.elements
                }

                if len(mol_ratio) > 1:
                    mix_enthalpy = self.tm.calc_mutinary_multilattice_mix_enthalpy(
                        mol_ratio=mol_ratio,
                        binary_dict=self.data,
                        end_member_dict=self.end_member,
                        correction=self.correction,
                        temperature=temperature,
                        transition_temperatures=self.transition_temperatures,
                        model="regular",
                    )
                    config_entropy = self.tm.calc_config_entropy(mol_ratio)

                    if isinstance(mix_enthalpy, dict):
                        entry.update_energy(
                            self.tm.calc_gibbs_energy(
                                mix_enthalpy[lattice], config_entropy, temperature
                            )
                            * entry.composition.num_atoms
                        )
                        new_pd_entries.append(entry)

        new_pd_entries = self.process_elements_hull(
            composition, temperature, new_pd_entries
        )
        return PhaseDiagram(new_pd_entries)

    def process_elements_hull(
        self,
        composition: List[str],
        temperature: float,
        pd_entries_list: List[PDEntryLocal],
    ) -> List[PDEntryLocal]:
        """
        Processes the elements in the composition to update the phase diagram entries.

        Args:
                composition (List[str]): A list of elements in the alloy system.
                temperature (float): The temperature in Kelvin.
                pd_entries_list (List[PDEntryLocal]): A list of PDEntryLocal objects representing the phase diagram entries.

        Returns:
                List[PDEntryLocal]: An updated list of PDEntryLocal objects.
        """
        for ele in composition:
            if ele in ["Fe", "Ti", "Mn", "Hf", "Zr"] and self.correction:
                for key, value in self.end_member[ele].items():
                    name = Composition(ele)
                    temp_energy = (
                        value
                        - value * temperature / self.transition_temperatures[ele][2]
                        if key == self.transition_temperatures[ele][1]
                        else value
                    )
                    pd_entries_list.append(
                        PDEntryLocal(
                            composition=name,
                            energy=temp_energy * name.num_atoms,
                            name=f"{name.alphabetical_formula}_{key}",
                        )
                    )
            else:
                for key, value in self.end_member[ele].items():
                    name = Composition(ele)
                    pd_entries_list.append(
                        PDEntryLocal(
                            composition=name,
                            energy=value * name.num_atoms,
                            name=f"{name.alphabetical_formula}_{key}",
                        )
                    )

        return pd_entries_list

    def find_enthalpy_entropy_composition(
        self,
        composition: List[str],
        lattice: str,
        mol_ratio: List[float],
        temperature: float,
    ) -> tuple[Union[int, float], float, dict[Any, Any]]:
        """
        Finds the enthalpy, entropy, and mole ratio for a given composition and lattice at a specific temperature.

        Args:
                composition (List[str]): A list of elements in the alloy.
                lattice (str): The lattice structure.
                mol_ratio (List[float]): A list of mole fractions.
                temperature (float): The temperature in Kelvin.

        Returns:
                Union[float, float, Dict[str, float]]: The enthalpy, entropy, and mole ratio for the alloy composition.
        """

        mol_ratio = dict(zip(composition, mol_ratio))
        mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}

        if len(mol_ratio) == 1:
            mix_enthalpy = 0
        else:
            mix_enthalpy = self.tm.calc_mutinary_multilattice_mix_enthalpy(
                mol_ratio=mol_ratio,
                binary_dict=self.data,
                end_member_dict=self.end_member,
                correction=self.correction,
                temperature=temperature,
                transition_temperatures=self.transition_temperatures,
                model="regular",
            )
            mix_enthalpy = (
                min(list(mix_enthalpy.values()))
                if lattice == "min"
                else mix_enthalpy[lattice]
            )

        entropy = self.tm.calc_config_entropy(mol_ratio)
        return mix_enthalpy, entropy, mol_ratio

    @staticmethod
    def make_PD_entry(mol_ratio: Dict[str, float], energy: float) -> PDEntryLocal:
        """
        Creates a phase diagram entry (PDEntryLocal) from the mole ratio and energy.

        Args:
                mol_ratio (Dict[str, float]): The mole ratio of elements.
                energy (float): The energy value.

        Returns:
                PDEntryLocal: A PDEntryLocal object representing a phase diagram entry.
        """
        name = Composition(Composition(mol_ratio).get_integer_formula_and_factor()[0])
        return PDEntryLocal(composition=name, energy=energy * name.num_atoms)
