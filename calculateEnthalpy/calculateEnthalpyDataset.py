import json
from mp_api.client import MPRester
from callMpAPI.utils import getAPIKey
from tqdm import tqdm
from calculateEnthalpy.helper_functions.thermo_math import thermo_maths
from calculateEnthalpy.helper_functions.grid_code import create_multinary, create_mol_grid
from calculateEnthalpy.helper_functions.data_utils import load_json, extract_ele_list
import os
import warnings


# warnings.filterwarnings("ignore")

class calculateEnthalpyDataset:

    def __init__(self,
                 folder_path: str = "./data/input_data",
                 lattice: str = "bcc",
                 source: str = "pravan",
                 max_alloy_n: int = 5,
                 im_flag: bool = False
                 ):

        #access the input folders
        self.source = source
        self.lattice = lattice
        self.im_flag = im_flag
        self.binary_dict = load_json(folder_path=folder_path, lattice=lattice, source=source)
        self.ele_list = extract_ele_list(folder_path=folder_path, lattice=lattice, source=source)
        self.max_alloy_n = max_alloy_n
        self.min_alloy_n = 2
        self.grid_size = 20
        self.tm = thermo_maths(binary_dict=self.binary_dict)
        self.output_folder_path = f"./data/output_data/"
        self.result_dict = {}
        self.mpr = MPRester(
            api_key=getAPIKey("/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/callMpAPI/api_key.txt"))
        self.fields = ['composition', 'formation_energy_per_atom', 'energy_above_hull', 'chemsys']

        self._create_output_folder()

    def _create_output_folder(self):
        lfolder = os.listdir(self.output_folder_path)
        count = 0
        for i in lfolder:
            if f"{self.source}_{self.lattice}" in i:
                count += 1
        self.out_file_path = f"{self.output_folder_path}/{self.source}_{self.lattice}_{count + 1}"
        os.mkdir(self.out_file_path)

    def _compute_multinaries(self):
        self.multinaries = create_multinary(element_list=self.ele_list,
                                            no_comb=list(range(self.min_alloy_n, self.max_alloy_n + 1)))

    def _extract_intermetallic(self,
                               comp):
        docs = self.mpr.materials.summary.search(chemsys=comp, fields=self.fields)
        # print(docs)
        return docs

    def intermetallic_sorter(self, n_nary_dict, docs):

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

    def find_n_nary_intermetallic(self, n_nary_dict):

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
                        # print(e)
                        continue  #except statement in case the binary is not in the dataset

                    entropy = self.tm.calc_configEntropy(mol_ratio=mol_ratio)
                    dim_alloy_dict[alloy]['-'.join([str(round(i, 2)) for i in mol_frac])] = {
                        'mix_enthalpy': mix_enthalpy,
                        'config_entropy': entropy
                    }

            if self.im_flag:
                self.result_dict[dim] = self.find_n_nary_intermetallic(n_nary_dict=dim_alloy_dict)

                with open(f'{self.out_file_path}/multinaries_im.json', 'w') as f:
                    json.dump(self.result_dict, f, ensure_ascii=False, indent=4)

            else:
                self.result_dict[dim] = dim_alloy_dict
                with open(f'{self.out_file_path}/multinaries_wo_im.json', 'w') as f:
                    json.dump(self.result_dict, f, ensure_ascii=False, indent=4)

        return 1


rc = calculateEnthalpyDataset(source="aziz",
                              lattice="bcc",
                              max_alloy_n=4).compute_mixing_enthalpies
