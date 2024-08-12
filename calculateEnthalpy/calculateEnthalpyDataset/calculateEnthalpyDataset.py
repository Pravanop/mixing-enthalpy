import json
import os
from typing import Union

from tqdm import tqdm

from mp_api.client import MPRester
from calculateEnthalpy.callMpAPI.utils import getAPIKey
from calculateEnthalpy.helper_functions.thermo_math import thermoMaths
from calculateEnthalpy.helper_functions.grid_code import create_multinary, create_mol_grid
from calculateEnthalpy.helper_functions.data_utils import DataUtils

import warnings
warnings.filterwarnings("ignore")

class calculateEnthalpyDataset:
    """
    Main class that implements the mixing enthalpy prediction model from DFT computed binary mixing enthalpies.
    This is an implementation of Zhang et al., A Fast and Robust Method for Predicting the Phase Stability of Refractory Complex
    Concentrated Alloys Using Pairwise Mixing Enthalpy http://dx.doi.org/10.2139/ssrn.4081906 with additional functionality added for off-equimolar compositions
    """
    def __init__(self,
                 folder_path: str,
                 lattice: str,
                 source: str,
                 max_alloy_n: int = 5,
                 im_flag: bool = False
                 ):
        """

        Args:
            folder_path: The folder path for the data folder. Recommended to create a Data/input data folder.
            lattice: lattice currently working on. Accepted values are 'bcc', 'fcc', 'hcp'
            source: a unique source identifier. The naming conventions for the input data can be read in the read me file for this folder.
            max_alloy_n: maximum alloy dimensionality to compute enthalpy for. Default is 5 (quinary)
            im_flag: whether intermetallics data should be added from Materials Project. Default is False.
        """

        #access the input folder
        self.source = source
        self.lattice = lattice
        self.im_flag = im_flag
        self.binary_dict = DataUtils().load_json(folder_path=folder_path, lattice=lattice, source=source)
        self.ele_list = DataUtils().extract_ele_list(folder_path=folder_path, lattice=lattice, source=source)

        self.max_alloy_n = max_alloy_n
        self.min_alloy_n = 2
        self.grid_size = 20

        self.tm = thermoMaths(binary_dict=self.binary_dict)
        self.output_folder_path = f"../../data/output_data/"
        self.result_dict = {}

        self.mpr = MPRester(
            api_key=getAPIKey("/callMpAPI/api_key.txt"))
        self.fields = ['composition', 'formation_energy_per_atom', 'energy_above_hull', 'chemsys']

        self._create_output_folder()

    def _create_output_folder(self) -> None:
        """
        Creates the output folder
        Returns: None

        """
        lfolder = os.listdir(self.output_folder_path)
        count = 0
        for i in lfolder:
            if f"{self.source}_{self.lattice}" in i:
                count += 1
        self.out_file_path = f"{self.output_folder_path}/{self.source}_{self.lattice}_{count + 1}"
        os.mkdir(self.out_file_path)

    def _compute_multinaries(self) -> None:
        """
        Computes all the multinaries possible ranging from 2 to self.max_alloy_n from the provided element list
        Returns: None

        """
        self.multinaries = create_multinary(element_list=self.ele_list,
                                            no_comb=list(range(self.min_alloy_n, self.max_alloy_n + 1)))

    def _extract_intermetallic(self,
                               comp: list[str]) -> Union[list, dict[list]]:
        """
        Extracts all the intermetallic properties for the given compositional space
        Args:
            comp: list of compositional elements of the form 'E1-E2-E3'

        Returns: A list of dicts containing information about the intermetallic properties

        """
        docs = self.mpr.materials.summary.search(chemsys=comp, fields=self.fields)
        return docs

    @staticmethod
    def intermetallic_sorter(
                             n_nary_dict : dict,
                             docs: Union[list, dict[list]]) -> dict:
        """
        Extracts all the intermetallic properties for the given compositional space and updates the dimensional dictionary
        Args:
            n_nary_dict: dictionary containing the composition and mixing enthalpy information
            docs: the search document from materials project

        Returns: the updated version of the dimensional dictionary with a key 'intermetallic' added to each entry

        """
        for idx, doc in enumerate(docs):

            temp_dict = {
                'formula_pretty': str(doc.composition.reduced_formula),
                'formation_energy_per_atom': doc.formation_energy_per_atom,
                'energy_above_hull': doc.energy_above_hull
            }
            comp = str(doc.chemsys)
            if "intermetallic" in n_nary_dict[comp]:
                n_nary_dict[comp]['intermetallic'].append(temp_dict)
            else:
                n_nary_dict[comp]['intermetallic'] = [temp_dict]

        return n_nary_dict

    def find_n_nary_intermetallic(self,
                                  n_nary_dict: dict) -> dict:
        """
        Extracts all the intermetallic properties for the given dimensional composition dictionary
        Args:
            n_nary_dict: dictionary containing the composition and mixing enthalpy information

        Returns: updated n_nary_dict with all the intermetallic properties with key 'intermetallic' added to each entry

        """

        comps = list(n_nary_dict.keys())
        docs = self._extract_intermetallic(
            comp=comps,
        )
        n_nary_dict = self.intermetallic_sorter(
            n_nary_dict=n_nary_dict,
            docs=docs
        )

        return n_nary_dict

    @property
    def compute_mixing_enthalpies(self) -> int:
        """
        Copmutes all the mixing enthalpy values in ev/atom and saves it in a json file that is stored in the specified output folder
        Returns: A return code of 1 implies successful computation

        """
        self._compute_multinaries()

        for dim, alloys in self.multinaries.items():
            mol_grid = create_mol_grid(dim, self.grid_size)
            dim_alloy_dict = {}
            for alloy_idx, alloy in enumerate(tqdm(alloys, desc=f"Accessing {dim}-nary")):
                ele_list = alloy.split('-')
                dim_alloy_dict[alloy] = {}
                for mol_idx, mol_frac in enumerate(mol_grid):
                    mol_ratio = dict(zip(ele_list, mol_frac))
                    mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}
                    try:
                        mix_enthalpy = self.tm.calc_multinary_mixEnthalpy(
                            alloy_comp=alloy,
                            mol_ratio=mol_ratio
                        )
                    except KeyError as e:
                        continue  #except statement in case the binary is not in the dataset, ideally the input data should be complete

                    entropy = self.tm.calc_configEntropy(mol_ratio=mol_ratio)
                    dim_alloy_dict[alloy]['-'.join([str(round(i, 2)) for i in mol_frac])] = {
                        'mix_enthalpy': mix_enthalpy,
                        'config_entropy': entropy
                    } #sorted strings for the alloy composition is important to mantain homogeneity.

            if self.im_flag:
                self.result_dict[dim] = self.find_n_nary_intermetallic(n_nary_dict=dim_alloy_dict)

                with open(f'{self.out_file_path}/multinaries_im.json', 'w') as f:
                    json.dump(self.result_dict, f, ensure_ascii=False, indent=4)

            else:
                self.result_dict[dim] = dim_alloy_dict
                with open(f'{self.out_file_path}/multinaries_wo_im.json', 'w') as f:
                    json.dump(self.result_dict, f, ensure_ascii=False, indent=4)

        return 1

