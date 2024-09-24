from typing import Union, Tuple, List, Dict, Iterable

import numpy as np
import pandas as pd
from numpy import ndarray
from pandas import DataFrame
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Element
from tqdm import tqdm

from far_heaa.data_processing.intermetallic_extraction import IntermetallicExtractions
from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.grids_and_combinations.grid_creation import CompositionGrid
from far_heaa.math_operations.thermo_calculations import ThermoMaths
from far_heaa.phase_diagram.convex_hull import ConvexHull


class GridIterator:

    def __init__(
        self,
        grid_size: int,
        tm: ThermoMaths,
        data: Dict,
        end_member: Dict,
        flags: Dict,
        api_key: str,
    ):
        """
        Initializes the class with the necessary parameters and sets up the ConvexHull instance.

        Args:
                grid_size (int): The size of the molar grid.
                tm (ThermoMaths): An instance of the ThermoMaths class for thermodynamic calculations.
                data (Dict): A dictionary containing binary mixing enthalpy data.
                end_member (Dict): A dictionary containing end member energies for various elements.
                flags (Dict): A dictionary of flags controlling behavior (e.g., 'im_flag', 'equi_flag', 'correction').
                api_key (str): The API key for accessing external data sources.
        """

        self.grid_size = grid_size
        self.end_member = end_member
        self.im_flag = flags["im_flag"]
        self.equi_flag = flags["equi_flag"]
        self.correction = flags["correction"]
        self.tm = tm
        self.api_key = api_key
        self.data = data
        self.convex_hull = ConvexHull(
            end_member=end_member,
            grid_size=grid_size,
            flags=flags,
            tm=tm,
            api_key=api_key,
            data=data,
        )

        self.temp_grid = None

    def upper_limit(
        self, composition: List[str], mol_ratio: List[float], phase_flag: bool
    ) -> None:
        """
        Calculates the upper temperature limit for the grid based on the average melting temperature.

        Args:
                composition (List[str]): A list of elements in the alloy.
                mol_ratio (List[float]): The mole fractions of the elements.
                phase_flag (bool): A flag indicating whether phase equilibrium should be considered.

        Returns:
                None: Updates the temperature grid (`self.temp_grid`) based on the phase flag and average melting temperature.
        """
        avg_tm = self.tm.avg_T_melt(composition, mol_ratio)
        if phase_flag:
            self.temp_grid = np.arange(0, avg_tm + 200, 20, dtype=np.float64)
        else:
            self.temp_grid = np.arange(0, avg_tm + 100, 100, dtype=np.float64)

    def temp_iterator(
        self, composition: List[str], temp_grid: List[float]
    ) -> Dict[float, PhaseDiagram]:
        """
        Iterates over a temperature grid and calculates the convex hull for each temperature.

        Args:
                composition (List[str]): A list of elements in the alloy system.
                temp_grid (List[float]): A list of temperatures to iterate over.

        Returns:
                Dict[float, PhaseDiagram]: A dictionary mapping each temperature to its corresponding PhaseDiagram.

        Example::

                results = self.temp_iterator(composition=['Fe', 'Ni'], temp_grid=[300, 500, 700])
                # Returns a dictionary with PhaseDiagram objects for each temperature.
        """
        n_alloy = len(composition)
        all_combs = MultinaryCombinations.create_multinary(
            element_list=composition, no_comb=list(range(2, n_alloy + 1))
        )

        im_list = []
        if self.im_flag:
            for dimensionality, alloy_list in all_combs.items():
                im_list += IntermetallicExtractions.get_MP_intermetallic(
                    alloy_list, api_key=self.api_key
                )

        PD_temp_comp_dict = {}
        for temp in tqdm(temp_grid, desc="Running Temperature"):
            PD_temp_comp_dict[temp] = self.convex_hull.make_convex_hull(
                temperature=temp, composition=composition, batch_tag=True, im=im_list
            )

        return PD_temp_comp_dict

    def misc_temperature_across_grid(
        self,
        composition: List[str],
        mol_grid_size: Union[int, ndarray],
        lattice: str,
        phase_flag: bool = False,
        **kwargs
    ) -> Tuple[ndarray, ndarray]:
        """

        Computes the miscibility temperature across a mole fraction grid for a given composition.

        Args:
                composition (List[str]): A list of elements in the alloy.
                mol_grid_size (Union[int, ndarray]): The size of the mole fraction grid. If not an integer, it is assumed to be the grid itself.
                lattice (str): The lattice structure to consider (e.g., 'FCC', 'BCC').
                phase_flag (bool, optional): If True, use finer temperature increments. Default is False.
                **kwargs: Additional arguments such as 'conv_hull' and 'temp_grid' if precomputed values are provided.

        Returns:
                Tuple[ndarray, ndarray]: A tuple containing:
                        - mol_grid (ndarray): The mole fraction grid for the alloy.
                        - stables (ndarray): The miscibility temperature for each composition in the grid.

        Example::
                mol_grid, stable_temps = self.misc_temperature_across_grid(
                        composition=['Fe', 'Ni'],
                        mol_grid_size=20,
                        lattice='FCC',
                        phase_flag=True
                )
                # Returns the mole fraction grid and miscibility temperatures across the grid.
        """

        n = len(composition)

        # Create the mole fraction grid
        if isinstance(mol_grid_size, int):
            mol_grid = CompositionGrid.create_mol_grid(n, mol_grid_size)
        else:
            mol_grid = mol_grid_size

        # Compute average melting temperature for each element in the composition
        temp_ele = []
        for i in composition:
            temp_ele.append(self.tm.avg_T_melt(i, list(mol_grid)))

        # Check if precomputed convex hull data is provided
        if "conv_hull" in kwargs:
            conv_hull_dict = kwargs["conv_hull"]
            temp_grid = kwargs["temp_grid"]
        else:
            # Generate the temperature grid based on the phase_flag
            if phase_flag:
                temp_grid = list(
                    np.arange(0, max(temp_ele) + 100, 20, dtype=np.float64)
                )
            else:
                temp_grid = list(
                    np.arange(0, max(temp_ele) + 100, 100, dtype=np.float64)
                )

            # Generate the convex hull at each temperature
            conv_hull_dict = self.temp_iterator(
                composition=composition, temp_grid=temp_grid
            )

        # Initialize stables array to store miscibility temperature
        stables = np.ones_like(mol_grid[:, 0]) * -1

        done_idx = []
        for idx_temp, temp in enumerate(temp_grid):

            for idx, mol in enumerate(mol_grid):
                if idx in done_idx:
                    continue

                # Calculate enthalpy and entropy for the current mole ratio
                mix_enthalpy, entropy, mol_ratio = (
                    self.convex_hull.find_enthalpy_entropy_composition(
                        mol_ratio=list(mol),
                        composition=composition,
                        temperature=temp,
                        lattice=lattice,
                    )
                )

                # Check if the composition is stable at this temperature
                is_stable = self.get_decomp_and_e_hull(
                    mol_ratio=mol_ratio,
                    temp=temp,
                    conv_hull=conv_hull_dict[temp],
                    mix_enthalpy=mix_enthalpy,
                    entropy=entropy,
                )

                if is_stable is not None:
                    # Check if the energy above the hull is within a tolerance of 1 meV/atom
                    if np.isclose(is_stable[1], 0.0, atol=1e-3):
                        stables[idx] = temp
                        done_idx.append(idx)
                    else:
                        continue
                else:
                    stables[idx] = temp
                    done_idx.append(idx)

        return mol_grid, stables

    def e_hull_across_grid(
        self,
        composition: List[str],
        mol_grid_size: Union[int, Iterable],
        lattice: str,
        single_temp_flag: bool = False,
        temp_gradation: int = 300,
    ) -> Tuple[ndarray, ndarray, ndarray]:
        """

        Computes the energy above the convex hull (E_hull) across a mole fraction grid for a given alloy composition.

        Args:
                composition (List[str]): A list of elements in the alloy.
                mol_grid_size (Union[int, Iterable]): The size of the mole fraction grid. If not an integer, it is assumed to be the grid itself.
                lattice (str): The lattice structure to consider (e.g., 'FCC', 'BCC').
                single_temp_flag (bool, optional): If True, computes the energy above the hull at a single temperature. Default is False.
                temp_gradation (int, optional): The temperature step size for the grid, in Kelvin. Default is 300 K.

        Returns:
                Tuple[ndarray, ndarray, ndarray]: A tuple containing:
                        - mol_grid (ndarray): The mole fraction grid for the alloy.
                        - stables (ndarray): The energy above the convex hull for each composition and temperature.
                        - temp_grid (ndarray): The grid of temperatures considered.

        Example::

                mol_grid, e_hull_values, temp_grid = self.e_hull_across_grid(
                        composition=['Fe', 'Ni'],
                        mol_grid_size=20,
                        lattice='FCC',
                        single_temp_flag=False,
                        temp_gradation=200
                )
                # Returns the mole fraction grid, energy above the hull values, and the temperature grid.
        """

        n = len(composition)

        # Create the mole fraction grid
        if isinstance(mol_grid_size, int):
            mol_grid = CompositionGrid.create_mol_grid(n, mol_grid_size)
        else:
            mol_grid = mol_grid_size

        # Compute average melting temperature for each element in the composition
        temp_ele = []
        for i in composition:
            temp_ele.append(self.tm.avg_T_melt(i, list(mol_grid)))

        # Create the temperature grid
        if single_temp_flag:
            temp_grid = [temp_gradation]
        else:
            temp_grid = list(
                np.arange(200, max(temp_ele) + 100, temp_gradation, dtype=np.float64)
            )

        # Generate the convex hull at each temperature
        conv_hull_dict = self.temp_iterator(
            composition=composition, temp_grid=temp_grid
        )

        # Initialize stables array to store energy above the convex hull
        stables = np.ones((len(mol_grid), len(temp_grid)))

        # Iterate over the temperature grid and mole fraction grid
        for idx_temp, temp in enumerate(temp_grid):
            for idx, mol in enumerate(mol_grid):

                # Calculate enthalpy and entropy for the current mole ratio
                mix_enthalpy, entropy, mol_ratio = (
                    self.convex_hull.find_enthalpy_entropy_composition(
                        mol_ratio=list(mol),
                        composition=composition,
                        temperature=temp,
                        lattice=lattice,
                    )
                )

                # Check if the composition is stable at this temperature
                is_stable = self.get_decomp_and_e_hull(
                    mol_ratio=mol_ratio,
                    temp=temp,
                    conv_hull=conv_hull_dict[temp],
                    mix_enthalpy=mix_enthalpy,
                    entropy=entropy,
                )

                # Store the energy above the convex hull or set it to 0.0 if stable
                if is_stable is not None:
                    if np.isclose(is_stable[1], 0.0, atol=1e-3):
                        stables[idx][idx_temp] = 0.0
                    else:
                        stables[idx][idx_temp] = is_stable[1]
                else:
                    stables[idx][idx_temp] = 0.0

        return mol_grid, stables, temp_grid

    def uni_molar_misc_temperature(
        self,
        mol_ratio: List[float],
        composition: List[str],
        lattice: str,
        phase_flag: bool = False,
        batch_tag: bool = False,
        **kwargs
    ) -> Union[float, Tuple[float, int]]:
        """
        Calculates the miscibility temperature for a given alloy composition and mole fraction.

        Args:
                mol_ratio (List[float]): The mole fractions of the elements in the alloy.
                composition (List[str]): A list of elements in the alloy.
                lattice (str): The lattice structure to consider (e.g., 'FCC', 'BCC').
                phase_flag (bool, optional): If True, applies a finer temperature increment. Default is False.
                batch_tag (bool, optional): If True, uses batch processing of intermetallic phases. Default is False.
                **kwargs: Additional arguments such as 'conv_hull' and 'im' (intermetallic data) if precomputed values are provided.

        Returns:
                Union[float, Tuple[float, int]]: The miscibility temperature for the given composition.
                        If no stable miscibility temperature is found, it returns a tuple (melting point, 0).

        Example:
                misc_temp = self.uni_molar_misc_temperature(
                        mol_ratio=[0.5, 0.5],
                        composition=['Fe', 'Ni'],
                        lattice='FCC',
                        phase_flag=True
                )
                # Returns the miscibility temperature for the alloy 'Fe-Ni' with a 50-50 mole ratio.
        """
        # Calculate the upper temperature limit based on the phase flag
        self.upper_limit(composition, mol_ratio, phase_flag)

        # Initialize intermetallic data based on batch processing flag
        if not batch_tag:
            n_alloy = len(composition)
            all_combs = MultinaryCombinations.create_multinary(
                element_list=composition, no_comb=list(range(2, n_alloy + 1))
            )
            im_list = []
            for dimensionality, alloy_list in all_combs.items():
                if self.im_flag:
                    im_list += IntermetallicExtractions.get_MP_intermetallic(
                        alloy_list, api_key=self.api_key
                    )
        else:
            im_list = kwargs.get("im", [])

        # Generate or retrieve the convex hull at 0K
        if "conv_hull" not in kwargs:
            conv_hull0 = self.convex_hull.make_convex_hull(
                temperature=0, composition=composition, batch_tag=True, im=im_list
            )
        else:
            conv_hull0 = kwargs["conv_hull"]

        # Iterate over the temperature grid
        for idx, temperature in enumerate(self.temp_grid):
            # Calculate enthalpy and entropy for the current mole ratio and temperature
            mix_enthalpy, entropy, mol_ratio_dict = (
                self.convex_hull.find_enthalpy_entropy_composition(
                    mol_ratio=mol_ratio,
                    composition=composition,
                    temperature=float(temperature),
                    lattice=lattice,
                )
            )

            # Update the convex hull for the current temperature
            conv_hull = self.convex_hull.update_convex_hull(
                composition=composition, conv_hull=conv_hull0, temperature=temperature
            )

            # Check if the composition is stable at this temperature
            is_stable = self.get_decomp_and_e_hull(
                mol_ratio=mol_ratio_dict,
                temp=float(temperature),
                conv_hull=conv_hull,
                mix_enthalpy=mix_enthalpy,
                entropy=entropy,
            )

            # Return the temperature if the alloy is stable at this temperature
            if is_stable is not None:
                if np.isclose(is_stable[1], 0.0, atol=1e-3):
                    return float(temperature)
                else:
                    continue
            else:
                return float(temperature)

        # If no stable temperature is found, return the upper temperature limit (melting point)
        return float(self.temp_grid[-1]), 0

    def find_decomp_products(
        self,
        mol_ratio: List[float],
        composition: List[str],
        temperature: float,
        lattice: str,
        batch_tag: bool = False,
        **kwargs
    ) -> Union[tuple[DataFrame, float], tuple[None, int]]:
        """

        Finds the decomposition products and energy above the convex hull for a given alloy composition
        and mole ratio at a specific temperature.

        Args:
                mol_ratio (List[float]): The mole fractions of the elements in the alloy.
                composition (List[str]): A list of elements in the alloy.
                temperature (float): The temperature in Kelvin for the calculations.
                lattice (str): The lattice structure to consider (e.g., 'FCC', 'BCC').
                batch_tag (bool, optional): If True, uses batch processing of intermetallic phases. Default is False.
                **kwargs: Additional arguments such as 'im' (intermetallic data) if precomputed values are provided.

        Returns:
                Union[Tuple[pd.DataFrame, float], int]: A tuple containing:
                        - pd.DataFrame: A DataFrame with the decomposition products and their fractions.
                        - float: The energy above the convex hull.
                        If no decomposition is found, it returns None and 0.

        Example:
                products, energy_above_hull = self.find_decomp_products(
                        mol_ratio=[0.5, 0.5],
                        composition=['Fe', 'Ni'],
                        temperature=1200.0,
                        lattice='FCC',
                        batch_tag=True,
                        im=precomputed_im_data
                )
                # Returns a DataFrame with decomposition products and the energy above the hull or None if there is no decomposition products
        """

        full_comp = composition

        # Find enthalpy, entropy, and mole ratio for the current configuration
        mix_enthalpy, entropy, mol_ratio = (
            self.convex_hull.find_enthalpy_entropy_composition(
                mol_ratio=mol_ratio,
                composition=composition,
                temperature=temperature,
                lattice=lattice,
            )
        )

        # Get intermetallic data if batch processing is enabled
        if batch_tag and self.im_flag:
            im_list = kwargs.get("im", [])
        else:
            n_alloy = len(composition)
            all_combs = MultinaryCombinations.create_multinary(
                element_list=composition, no_comb=list(range(2, n_alloy + 1))
            )
            im_list = []
            for dimensionality, alloy_list in all_combs.items():
                if self.im_flag:
                    im_list += IntermetallicExtractions.get_MP_intermetallic(
                        alloy_list, api_key=self.api_key
                    )

        # Calculate configurational entropy
        entropy = self.tm.calc_config_entropy(mol_ratio)
        composition = list(mol_ratio.keys())

        # Exclude certain intermetallic phases based on the full composition
        exclude = list(set(full_comp).difference(set(composition)))
        im_list_exclude = []
        if exclude:
            for j in exclude:
                for idx, i in enumerate(im_list):
                    if Element(j) in i.elements:
                        im_list_exclude.append(i)

        im_list = list(set(im_list).difference(im_list_exclude))

        # Create the convex hull for the current composition and temperature
        if batch_tag:
            conv_hull = self.convex_hull.make_convex_hull(
                temperature=float(temperature),
                composition=composition,
                batch_tag=True,
                im=im_list,
            )
        else:
            conv_hull = self.convex_hull.make_convex_hull(
                temperature=float(temperature), composition=composition
            )

        # Check if the composition is stable at this temperature
        is_stable = self.get_decomp_and_e_hull(
            mol_ratio=mol_ratio,
            temp=float(temperature),
            conv_hull=conv_hull,
            mix_enthalpy=mix_enthalpy,
            entropy=entropy,
        )

        # If stable, return the decomposition products and energy above the convex hull
        if is_stable:
            col1, col2 = [], []
            for key, value in is_stable[0].items():
                col1.append(key.name)
                col2.append(round(value, 3))

            df = pd.DataFrame({"Product": col1, "Fraction": col2})
            return df, np.round(is_stable[1], 4)

        # If no decomposition is found, return None and 0
        return None, 0

    def get_decomp_and_e_hull(
        self,
        mol_ratio: Dict[str, float],
        temp: float,
        conv_hull: PhaseDiagram,
        mix_enthalpy: float,
        entropy: float,
    ) -> Union[Tuple[Dict, float], None]:
        """
        Calculates the Gibbs free energy and determines the decomposition products and energy above the convex hull.

        Args:
                mol_ratio (Dict[str, float]): A dictionary of mole ratios for the alloy.
                temp (float): The temperature in Kelvin.
                conv_hull (PhaseDiagram): The convex hull for the given temperature and composition.
                mix_enthalpy (float): The mixing enthalpy for the composition.
                entropy (float): The configurational entropy for the composition.

        Returns:
                Union[Tuple[Dict, float], None]: A tuple containing:
                        - Dict: Decomposition products and their mole fractions.
                        - float: Energy above the convex hull.
                        Returns None if no decomposition is found.

        Example:
                decomp_products, energy_above_hull = self.get_decomp_and_e_hull(
                        mol_ratio={'Fe': 0.5, 'Ni': 0.5},
                        temp=1200.0,
                        conv_hull=convex_hull,
                        mix_enthalpy=-0.5,
                        entropy=0.1
                )
                # Returns the decomposition products and energy above the convex hull.
        """

        # Calculate Gibbs free energy at the given temperature
        gibbs = self.tm.calc_gibbs_energy(
            enthalpy=mix_enthalpy, entropy=entropy, temperature=temp
        )

        # Create a phase diagram entry for the given mole ratio and Gibbs free energy
        pdEntry = self.convex_hull.make_PD_entry(energy=gibbs, mol_ratio=mol_ratio)

        try:
            # Get the decomposition products and energy above the convex hull
            answer = conv_hull.get_decomp_and_e_above_hull(
                entry=pdEntry, check_stable=False
            )
            return answer

        except ValueError:
            # Return None if no decomposition is found
            return None
