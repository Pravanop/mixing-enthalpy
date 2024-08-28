import itertools
import pickle
from typing import Tuple
import numpy as np
from tqdm import tqdm

from calculateEnthalpy.helper_functions.data_utils import DataUtils
from calculateEnthalpy.helper_functions.grid_code import create_multinary
from calculateEnthalpy.helper_functions.thermo_math import thermoMaths
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram


class reactionPathways:
    """
    Calculates deposition sequence for a certain quartenary HEA based on mixing parameter governed by predicted from miscibility temperature
    """
    def __init__(self,
                 lattice: str,
                 source: str,
                 inp_meta: dict)-> None:
        """

        Args:
            lattice:  lattice currently working on. Accepted values are 'bcc', 'fcc', 'hcp'
            source: a unique source identifier. The naming conventions for the input data can be read in the read me file for this folder
            inp_meta: A dictionary containing melting temperature information of the base elements. An example is below
        """
        self.lattice = lattice
        self.source = source
        self.inp_meta = inp_meta

        self.ele_list_main = list(inp_dict.keys())

        binary_dict = DataUtils.load_json(folder_path=f"../../data/input_data/", lattice=lattice, source=source)
        self.tm = thermoMaths(binary_dict=binary_dict)
        self.pD = phaseDiagram(source=source, lattice=lattice, version_no=1, im_flag=False, abs_file_path=f"../../data/output_data", binary_dict_path=f"../../data/input_data")

        self.binary_mol = [[0.8, 0.2], [0.7, 0.3], [0.5, 0.5]]
        self.ternary_mol = [[0.45, 0.45, 0.1], [0.4, 0.4, 0.2], [0.33, 0.33, 0.34]]
        self.misc_T_scores = None
        self.enthalpy_scores = None

        self._initialize_for_phaseSpace()

    def _update_pathways_dict(self,
                              pathways_dict: dict)-> None:
        """

        Args:
            pathways_dict: Sets a pathways dict incase one has been precomputed and stored elsewhere.

        Returns:

        """

        self.pathway_energies_temp = pathways_dict

    def _compute_all_comps(self)-> dict:
        """
        Creates all subset HEAs within the given composition.
        Returns: a dictionary of the form
	{
	3: ['E1-E2-E3', 'E1-E3-E4'....],
	2: ....
	}

        """
        return create_multinary(self.ele_list_main, no_comb=list(range(2, len(self.ele_list_main) + 1)))

    def _compute_all_pathways(self) -> list:
        """
        Computes all the deposition sequences possible for the given composition
        Returns: List of permutations in the form [(1, 2, 3,4), (1, 3, 2, 4), (2, 3, 1, 4)....]

        """
        return list(itertools.permutations(self.ele_list_main, len(self.ele_list_main)))

    def _compute_x_coord(self) -> None:
        """
        Computes the x coordinates for the subset alloys. Assigns increasing x value for higher dim alloys.
        Returns: None

        """
        coords_dict = {}
        for comp, _ in self.enthalpy_dict.items():
            if '-' in comp:
                length = len(comp.split('-'))
            else:
                length = 0

            coords_dict[comp] = length * 1.5
            for i in comp.split('-'):
                coords_dict[comp] += self.inp_meta[i]['coords'] # Base coords are present in the input meta dictionary

        self.coords_dict = coords_dict

    def _compute_enthalpies_and_temperatures_for_each_sub_alloy(self) -> None:
        """
        The main function that calculates the mixing enthalpies, and gets the miscibility temperature from the created convex hull.
        Returns: None

        """

        enthalpy_dict = {}
        misc_temp_dict = {}
        sub_enthalpy_dict = {}
        sub_misc_temp_dict = {}
        for key, value in tqdm(self.all_comps.items(), desc="Calculating T_m for phase space"):
            # iterates through the various dims
            for idx, alloy in enumerate(value):
                # iterates through the alloys in a dimensionality
                subset_ele_list = alloy.split("-")
                # hardcoded off equimolar ratios for each subset. Can be changed to better code later.
                if len(subset_ele_list) == 2:
                    for idx, mol in enumerate(self.binary_mol):

                        mol_ratio = mol
                        misc_temp = self.pD.find_misc_temperature(mol_ratio=mol_ratio,
                                                                  composition=alloy)
                        mol_ratio = dict(zip(subset_ele_list, mol_ratio))
                        alloy_key = '-'.join([f"{key},{value}" for key, value in mol_ratio.items()])
                        # print(misc_temp, alloy_key)
                        if idx == len(self.binary_mol) - 1:
                            enthalpy_dict[alloy] = self.tm.calc_multinary_mixEnthalpy(alloy_comp=alloy,
                                                                                      mol_ratio=mol_ratio)
                            if misc_temp is not None:
                                misc_temp_dict[alloy] = misc_temp
                            else:
                                misc_temp_dict[alloy] = 5000 #a very high temperature so that sorting, and subsequent functions don't get an error.
                        else:
                            # the offequimolar ratios are kept in a separate function to facilitate easier processing for plotting and visualization
                            sub_enthalpy_dict[alloy_key] = self.tm.calc_multinary_mixEnthalpy(alloy_comp=alloy,
                                                                                              mol_ratio=mol_ratio)

                            if misc_temp is not None:
                                sub_misc_temp_dict[alloy_key] = misc_temp
                            else:
                                sub_misc_temp_dict[alloy_key] = 5000

                elif len(subset_ele_list) == 3:
                    for idx, mol in enumerate(self.ternary_mol):
                        mol_ratio = mol
                        misc_temp = self.pD.find_misc_temperature(mol_ratio=mol_ratio,
                                                                  composition=alloy)
                        mol_ratio = dict(zip(subset_ele_list, mol_ratio))
                        alloy_key = '-'.join([f"{key},{value}" for key, value in mol_ratio.items()])
                        # print(misc_temp, alloy_key)
                        if idx == len(self.ternary_mol) - 1:

                            enthalpy_dict[alloy] = self.tm.calc_multinary_mixEnthalpy(alloy_comp=alloy,
                                                                                      mol_ratio=mol_ratio)
                            if misc_temp is not None:
                                misc_temp_dict[alloy] = misc_temp
                            else:
                                misc_temp_dict[alloy] = 5000
                        else:
                            sub_enthalpy_dict[alloy_key] = self.tm.calc_multinary_mixEnthalpy(alloy_comp=alloy,
                                                                                              mol_ratio=mol_ratio)
                            if misc_temp is not None:
                                sub_misc_temp_dict[alloy_key] = misc_temp
                            else:
                                sub_misc_temp_dict[alloy_key] = 5000

                else:
                    #hardcoded for quarternary
                    mol_ratio = [1 / len(subset_ele_list)] * len(subset_ele_list)
                    misc_temp = self.pD.find_misc_temperature(mol_ratio=mol_ratio,
                                                              composition=alloy)
                    mol_ratio = dict(zip(subset_ele_list, mol_ratio))
                    enthalpy_dict[alloy] = self.tm.calc_multinary_mixEnthalpy(alloy_comp=alloy,
                                                                              mol_ratio=mol_ratio)

                    if misc_temp is not None:
                        misc_temp_dict[alloy] = misc_temp
                    else:
                        misc_temp_dict[alloy] = 5000

        for key, value in self.inp_meta.items():
            misc_temp_dict.update({
                key: value['melting_temp']
            })
            enthalpy_dict.update({
                key: 0
            })

        self.enthalpy_dict = enthalpy_dict
        self.misc_temp_dict = misc_temp_dict
        self.sub_enthalpy_dict = sub_enthalpy_dict
        self.sub_misc_temp_dict = sub_misc_temp_dict

    def _initialize_for_phaseSpace(self) -> None:
        """
        Wrapper function that computes all the pre-processing code.
        Returns: None

        """
        self.all_comps = self._compute_all_comps()
        self.all_pathways = self._compute_all_pathways()
        self._compute_enthalpies_and_temperatures_for_each_sub_alloy()
        self._compute_x_coord()

    @staticmethod
    def compute_pathway_score(pathway: list) -> float:
        """
        The pathway score is the average miscible temperature for that path.
        Args:
            pathway: list of miscibility temperatures

        Returns: Average miscible temperature.

        """
        return sum(pathway) / len(pathway)

    def rank_pathway(self,
                     pathway_energies_temp: dict) -> Tuple[np.ndarray, np.ndarray]:
        """
        Ranks the pathways according to their average miscibility temperature.
        Args:
            pathway_energies_temp: The dictionary computed from compute_pathway_energies_temp

        Returns: The argsort for the enthalpy and miscibility temperatures to give rank later on

        """
        enthalpy_scores = []
        misc_T_scores = []
        for paths, values in pathway_energies_temp.items():
            temperatures = [i['misc_T'] for i in list(values.values())]
            enthalpies = [i['enthalpy'] for i in list(values.values())]
            sub_temperatures = [list(i['subs'].values()) for i in list(values.values()) if i['subs'] != {}]
            sub_temperatures = list(np.array(sub_temperatures).reshape(-1, 1))
            temperatures += sub_temperatures
            enthalpy_score = self.compute_pathway_score(pathway=enthalpies)
            misc_T_score = self.compute_pathway_score(pathway=temperatures)
            enthalpy_scores.append(enthalpy_score)
            misc_T_scores.append(misc_T_score)

        return np.argsort(enthalpy_scores), np.argsort(misc_T_scores)

    def _sub_match_key(self,
                       keys: list,
                       match_key: str) -> list:
        """
        Helper function to find the indices of the sub nodes in the node dictionary, to match off equimolar compositions with their equimolar counterparts.
        Args:
            keys: List of equimolar keys
            match_key: off-equimolar "sub-key"

        Returns: list of matching indices

        """
        eles = match_key.split('-')
        indices = []
        for idx, i in enumerate(keys):
            flag = 0
            for j in eles:
                if j in i:
                    flag += 1
                else:
                    continue

            if flag == len(eles) and len(eles) == len(i.split('-')):
                indices.append(idx)

        return indices

    @property
    def compute_pathway_energies_temp(self) -> dict:
        """
        Computes all the pathways with their node energies, and misc temperatures, along with sub-node information
        Returns: dictionary of pathway energies and misc temperatures in the form
        {'pathway1': {
                        alloy_name: {
                                    enthalpy:,
                                    misc_T:,
                                    sub: {
                                            sub_alloy1 : misc_T,
                                            sub_alloy2: ...,
                                         }
                                    }
                        },
        'pathway2': ...}

        """
        pathway_energies_temp = {}
        sub_keys = list(self.sub_misc_temp_dict.keys())
        for paths in self.all_pathways:
            alloys = []
            for i in range(len(self.ele_list_main)):
                if i == 0:
                    alloys.append(paths[0])
                else:
                    alloys.append('-'.join(paths[:i + 1]))
            pathway_energies_temp['-'.join(paths)] = {}
            for alloy in alloys:
                sorted_alloy = '-'.join(sorted(alloy.split('-')))
                indices = self._sub_match_key(keys=sub_keys, match_key=sorted_alloy)
                sub_dict = {}
                if indices:
                    for i in indices:
                        sub_dict[sub_keys[i]] = self.sub_misc_temp_dict[sub_keys[i]]

                pathway_energies_temp['-'.join(paths)].update({
                    alloy: {
                        'enthalpy': self.enthalpy_dict[sorted_alloy],
                        'misc_T': self.misc_temp_dict[sorted_alloy],
                        'subs': sub_dict
                    },
                })

        return pathway_energies_temp

    @property
    def pickler(self) -> dict:
        """
        Makes a dictionary object that stores all the information needed for plotting and visulization
        Returns: dictionary of needed information needed for plotting and visualization by plotReactionPathways

        """
        return {
            'all_pathways': self.all_pathways,
            'ele_list_main': self.ele_list_main,
            'all_comps': self.all_comps,
            'enthalpy_dict': self.enthalpy_dict,
            'misc_temp_dict': self.misc_temp_dict,
            'sub_misc_temp_dict': self.sub_misc_temp_dict,
            'sub_enthalpy_dict': self.sub_enthalpy_dict,
        }


#Example Runs.

if __name__ == '__main__':
    inp_dict = {
        'Ti': {
            # 'melting_temp': 2750,
            'melting_temp': 1950,
            'mols': 1,
            'coords': 0.2,
            'color': '#004488'
        },
        'Cr': {
            # 'melting_temp': 2183,
            'melting_temp': 2750,
            'mols': 1,
            'coords': 0.4,
            'color': "#DDAA33"
        },
        'W': {
            'melting_temp': 2128,
            'mols': 1,
            'coords': 0.6,
            'color': "#BB5566"
        },
        'Ta': {
            'melting_temp': 3293,
            'mols': 1,
            'coords': 0.8,
            'color': "#EE7733"
        }
    }
    lattice = "bcc"
    source = "pravan"

    rP = reactionPathways(lattice=lattice, source=source, inp_meta=inp_dict)
    pathway_energies_temp = rP.compute_pathway_energies_temp
    pathway_scores = rP.rank_pathway(pathway_energies_temp)

    pickle.dump((
        rP.pickler,
        pathway_energies_temp,
        pathway_scores),
        open(f'../../data/output_data/{source}_{lattice}_1'
			 f'/pathway_energies_temp_{"-".join(inp_dict.keys())}.p', 'wb'))
