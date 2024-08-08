import itertools
import pickle
from typing import Tuple

import numpy as np
from calculateEnthalpy.helper_functions.data_utils import load_json
from calculateEnthalpy.helper_functions.grid_code import create_multinary
from calculateEnthalpy.helper_functions.thermo_math import thermo_maths
from calculateEnthalpy.helper_functions.phaseDiagram import phaseDiagram
import matplotlib.pyplot as plt

inp_dict = {
    'Nb': {
        'melting_temp': 2750,
        'mols': 1,
        'coords': 0.2,
        'color': '#004488'
    },
    'V': {
        'melting_temp': 2183,
        'mols': 1,
        'coords': 0.4,
        'color': "#DDAA33"
    },
    'Zr': {
        'melting_temp': 2128,
        'mols': 1,
        'coords': 0.6,
        'color': "#BB5566"
    },
    'Ti': {
        'melting_temp': 1941,
        'mols': 1,
        'coords': 0.8,
        'color': "#EE7733"
    }
}
lattice = "bcc"
source = "aziz"


class reactionPathways:

    def __init__(self, lattice, source, inp_meta):

        self.lattice = lattice
        self.source = source
        self.inp_meta = inp_meta

        self.ele_list_main = list(inp_dict.keys())

        binary_dict = load_json(folder_path=f"./data/input_data/", lattice=lattice, source=source)
        self.tm = thermo_maths(binary_dict=binary_dict)
        self.pD = phaseDiagram(source=source, lattice=lattice, version_no=2)

        self.misc_T_scores = None
        self.enthalpy_scores = None

        self._initialize_for_phaseSpace()

    def _update_pathways_dict(self, pathways_dict):
        self.pathway_energies_temp = pathways_dict

    def _compute_all_comps(self):
        return create_multinary(self.ele_list_main, no_comb=list(range(2, len(self.ele_list_main) + 1)))

    def _compute_all_pathways(self):
        return list(itertools.permutations(self.ele_list_main, len(self.ele_list_main)))

    def _compute_x_coord(self):
        coords_dict = {}
        for comp, _ in self.enthalpy_dict.items():
            if '-' in comp:
                length = len(comp.split('-'))
            else:
                length = 0

            coords_dict[comp] = length * 1.5
            for i in comp.split('-'):
                coords_dict[comp] += self.inp_meta[i]['coords']

        self.coords_dict = coords_dict

    def _compute_enthalpies_and_temperatures_for_each_sub_alloy(self):

        enthalpy_dict = {}
        misc_temp_dict = {}
        for key, value in self.all_comps.items():
            for idx, alloy in enumerate(value):
                subset_ele_list = alloy.split("-")
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

    def _initialize_for_phaseSpace(self):
        self.all_comps = self._compute_all_comps()
        self.all_pathways = self._compute_all_pathways()
        self._compute_enthalpies_and_temperatures_for_each_sub_alloy()
        self._compute_x_coord()

    @staticmethod
    def compute_pathway_score(pathway: list):
        return sum(pathway) / len(pathway)

    def rank_pathway(self, pathway_energies_temp) -> Tuple[np.ndarray, np.ndarray]:
        enthalpy_scores = []
        misc_T_scores = []
        for paths, values in pathway_energies_temp.items():
            temperatures = [i['misc_T'] for i in list(values.values())]
            enthalpies = [i['enthalpy'] for i in list(values.values())]
            enthalpy_score = self.compute_pathway_score(pathway=enthalpies)
            misc_T_score = self.compute_pathway_score(pathway=temperatures)
            enthalpy_scores.append(enthalpy_score)
            misc_T_scores.append(misc_T_score)

        return np.argsort(enthalpy_scores), np.argsort(misc_T_scores)

    @property
    def compute_pathway_energies_temp(self) -> dict:
        pathway_energies_temp = {}

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
                pathway_energies_temp['-'.join(paths)].update({
                    alloy: {
                        'enthalpy': self.enthalpy_dict[sorted_alloy],
                        'misc_T': self.misc_temp_dict[sorted_alloy]
                    }
                })

        return pathway_energies_temp

    @property
    def pickler(self) -> dict:
        return {
            'all_pathways': self.all_pathways,
            'ele_list_main': self.ele_list_main,
            'all_comps': self.all_comps,
            'enthalpy_dict': self.enthalpy_dict,
            'misc_temp_dict': self.misc_temp_dict,
        }


rP = reactionPathways(lattice=lattice, source=source, inp_meta=inp_dict)
pathway_energies_temp = rP.compute_pathway_energies_temp
pathway_scores = rP.rank_pathway(pathway_energies_temp)

pickle.dump((
    rP.pickler,
    pathway_energies_temp,
    pathway_scores),
    open('/calculateEnthalpy/data/output_data/aziz_bcc_2/pathway_energies_temp.p', 'wb'))
