import json
import os
import time
from typing import Union, Tuple

import numpy as np
from tqdm import tqdm

from mp_api.client import MPRester
from calculateEnthalpy.callMpAPI.utils import getAPIKey
from calculateEnthalpy.helper_functions.thermo_math import thermoMaths
from calculateEnthalpy.helper_functions.grid_code import create_multinary, create_mol_grid
from calculateEnthalpy.helper_functions.data_utils import DataUtils
from emmet.core.thermo import ThermoType
import warnings
warnings.filterwarnings("ignore")

class calculateEnthalpyDataset:
    """
    Main class that implements the mixing enthalpy prediction model from DFT computed binary mixing enthalpies.
    This is an implementation of Zhang et al., A Fast and Robust Method for Predicting the Phase Stability of Refractory Complex
    Concentrated Alloys Using Pairwise Mixing Enthalpy http://dx.doi.org/10.2139/ssrn.4081906 with additional functionality added for off-equimolar compositions
    """
    def __init__(self,
                 input_folder_path: str,
                 lattice: str,
                 source: str,
                 max_alloy_n: int = 5,
                 im_flag: bool = False,
                 ele_list : Union[list, None] = None
                 ):
        """

        Args:
            input_folder_path: The folder path for the data folder. Recommended to create a Data/input data folder.
            lattice: lattice currently working on. Accepted values are 'bcc', 'fcc', 'hcp'
            source: a unique source identifier. The naming conventions for the input data can be read in the read me file for this folder.
            max_alloy_n: maximum alloy dimensionality to compute enthalpy for. Default is 5 (quinary)
            im_flag: whether intermetallics data should be added from Materials Project. Default is False.
        """

        #access the input folder
        self.source = source
        self.lattice = lattice
        self.im_flag = im_flag
        self.input_folder_path = input_folder_path
        # self.binary_dict = DataUtils().load_json(folder_path=folder_path, lattice=lattice, source=source)
        if ele_list is None:
            self.ele_list = DataUtils().extract_ele_list(folder_path=input_folder_path, lattice=lattice, source=source)
        else:
            self.ele_list = ele_list

        self.max_alloy_n = max_alloy_n
        self.min_alloy_n = 5
        self.grid_size = 20

        # self.tm = thermoMaths(binary_dict=self.binary_dict)
        self.output_folder_path = f"../../data/output_data/"
        self.result_dict = {}
        self.tm = thermoMaths()
        self.mpr = MPRester(
            api_key=getAPIKey("../callMpAPI/api_key.txt"))
        self.fields = ['composition', 'formation_energy_per_atom', 'energy_above_hull', 'chemsys']

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

    def process_binaries(self) -> Tuple[dict, list]:
        source_data_path = f"{self.input_folder_path}/{self.source}"
        lfolder = os.listdir(source_data_path)
        pairs = create_multinary(self.ele_list, no_comb=[2])
        master_binary = {}
        for file in lfolder:
            if file.endswith(".json"):
                lattice, source = file.split("_")
                master_binary[lattice] = DataUtils.load_json(folder_path=f"{self.input_folder_path}",
                                                  lattice=lattice,
                                                  source=self.source)

        result_dict = {}
        error_list = []
        lattice_value_pair = []
        for dimension, pair_list in pairs.items():
            for idx, pair in enumerate(tqdm(pair_list,  desc="Processing pairs created from element_list")):
                try:
                    lattice_value_pair = np.array([[key, master_binary[key][pair]] for key in master_binary.keys()])
                except KeyError:
                    error_list.append(pair)
                if len(lattice_value_pair) > 0:
                    result_dict[pair] = dict(zip(lattice_value_pair[:,0], np.round(lattice_value_pair[:,1].astype(float), 5)))

        return result_dict, error_list

    def store_outputs(self, results_dict) -> None:

        with open(f'{self.out_file_path}/all_lattices_binaries.json', 'w') as f:
            json.dump(results_dict, f, ensure_ascii=False, indent=4)

        print(f"All lattices binaries stored in {self.out_file_path}")

    @property
    def process_input_main(self):
        self._create_output_folder()
        results_dict, error_list = self.process_binaries()
        print("Certain pairs are not present in the input file: ", error_list)
        self.store_outputs(results_dict)
        return 1

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
            alloy_dict_values = []
            for alloy_idx, alloy in enumerate(tqdm(alloys, desc=f"Accessing {dim}-nary")):
                ele_list = alloy.split('-')
                alloy_mol_dict_values = []
                for mol_idx, mol_frac in enumerate(mol_grid):
                    mol_ratio = dict(zip(ele_list, mol_frac))
                    mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}
                    if len(mol_ratio) == 1:
                        continue

                    try:
                        # start = time.time()
                        mix_enthalpy = self.tm.calc_multinary_mixEnthalpy(
                            alloy_comp=alloy,
                            mol_ratio=mol_ratio
                        )
                        # stop = time.time()

                    except KeyError as e:
                        continue  #except statement in case the binary is not in the dataset, ideally the input data should be complete

                    # print(stop-start)
                    alloy_mol_dict_values.append(['-'.join([str(round(i, 2)) for i in mol_frac]),mix_enthalpy])

                if len(alloy_mol_dict_values) > 0:
                    alloy_dict_values.append([alloy, alloy_mol_dict_values])

            for j in alloy_dict_values:
                dim_alloy_dict[j[0]] = j[1]

            if self.im_flag:
                dim_alloy_dict = self.find_n_nary_intermetallic(n_nary_dict=dim_alloy_dict)

                with open(f'{self.out_file_path}/multinaries_im.json', 'w') as f:
                    json.dump(dim_alloy_dict, f, ensure_ascii=False, indent=4)

            else:
                with open(f'{self.out_file_path}/multinaries_wo_im_{dim}.json', 'w') as f:
                    json.dump(dim_alloy_dict, f, ensure_ascii=False, indent=4)

            # for n in range(self.min_alloy_n, self.max_alloy_n+1):
            #

        return 1


ced = calculateEnthalpyDataset(input_folder_path="../../data/input_data/",
                               lattice="bcc",
                               source="bokas2",
                               im_flag=True
                               # ele_list=['Cr', 'Fe','Mn','Ta','Ti','W','V','Mo','Nb','Zr']
                               )
print(ced.compute_mixing_enthalpies)