import itertools
import pickle
from typing import Tuple

import numpy as np
from calculateEnthalpy.helper_functions.data_utils import load_json
from calculateEnthalpy.helper_functions.grid_code import create_multinary, create_mol_grid
from calculateEnthalpy.helper_functions.thermo_math import thermo_maths
from calculateEnthalpy.helper_functions.phaseDiagram import phaseDiagram
import matplotlib.pyplot as plt


class reactionPathways:

    def __init__(self, lattice, source, inp_meta):

        self.lattice = lattice
        self.source = source
        self.inp_meta = inp_meta

        self.ele_list_main = list(inp_dict.keys())

        binary_dict = load_json(folder_path=f"../../data/input_data/", lattice=lattice, source=source)
        self.tm = thermo_maths(binary_dict=binary_dict)
        self.pD = phaseDiagram(source=source, lattice=lattice, version_no=2)

        self.binary_mol = [[0.8, 0.2], [0.7, 0.3], [0.5, 0.5]]
        self.ternary_mol = [[0.45, 0.45, 0.1], [0.4, 0.4, 0.2], [0.33, 0.33, 0.34]]
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
        sub_enthalpy_dict = {}
        sub_misc_temp_dict = {}
        print(self.all_comps)
        for key, value in self.all_comps.items():
            for idx, alloy in enumerate(value):
                print(alloy)
                subset_ele_list = alloy.split("-")
                if len(subset_ele_list) == 2:
                    for idx, mol in enumerate(self.binary_mol):

                        mol_ratio = mol
                        misc_temp = self.pD.find_misc_temperature(mol_ratio=mol_ratio,
                                                                  composition=alloy)
                        mol_ratio = dict(zip(subset_ele_list, mol_ratio))
                        alloy_key = '-'.join([f"{key},{value}" for key, value in mol_ratio.items()])
                        if idx == len(self.binary_mol) - 1:
                            enthalpy_dict[alloy] = self.tm.calc_multinary_mixEnthalpy(alloy_comp=alloy,
                                                                                      mol_ratio=mol_ratio)
                            print(enthalpy_dict[alloy])
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

                elif len(subset_ele_list) == 3:
                    for idx, mol in enumerate(self.ternary_mol):
                        mol_ratio = mol
                        misc_temp = self.pD.find_misc_temperature(mol_ratio=mol_ratio,
                                                                  composition=alloy)
                        mol_ratio = dict(zip(subset_ele_list, mol_ratio))
                        alloy_key = '-'.join([f"{key},{value}" for key, value in mol_ratio.items()])
                        if idx == len(self.ternary_mol) - 1:

                            enthalpy_dict[alloy] = self.tm.calc_multinary_mixEnthalpy(alloy_comp=alloy,
                                                                                      mol_ratio=mol_ratio)
                            print(enthalpy_dict[alloy])
                            if misc_temp is not None:
                                misc_temp_dict[alloy] = misc_temp
                            else:
                                misc_temp_dict[alloy] = 5000
                        else:
                            sub_enthalpy_dict[alloy_key] = self.tm.calc_multinary_mixEnthalpy(alloy_comp=alloy,
                                                                                              mol_ratio=mol_ratio)
                            print(sub_enthalpy_dict[alloy_key])
                            if misc_temp is not None:
                                sub_misc_temp_dict[alloy_key] = misc_temp
                            else:
                                sub_misc_temp_dict[alloy_key] = 5000

                else:
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
            sub_temperatures = [list(i['subs'].values()) for i in list(values.values()) if i['subs'] != {}]
            sub_temperatures = list(np.array(sub_temperatures).reshape(-1, 1))
            temperatures += sub_temperatures
            enthalpy_score = self.compute_pathway_score(pathway=enthalpies)
            misc_T_score = self.compute_pathway_score(pathway=temperatures)
            enthalpy_scores.append(enthalpy_score)
            misc_T_scores.append(misc_T_score)

        return np.argsort(enthalpy_scores), np.argsort(misc_T_scores)

    def _sub_match_key(self, keys, match_key):
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
        return {
            'all_pathways': self.all_pathways,
            'ele_list_main': self.ele_list_main,
            'all_comps': self.all_comps,
            'enthalpy_dict': self.enthalpy_dict,
            'misc_temp_dict': self.misc_temp_dict,
            'sub_misc_temp_dict': self.sub_misc_temp_dict,
            'sub_enthalpy_dict': self.sub_enthalpy_dict,
        }


if __name__ == '__main__':
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

    rP = reactionPathways(lattice=lattice, source=source, inp_meta=inp_dict)
    pathway_energies_temp = rP.compute_pathway_energies_temp
    pathway_scores = rP.rank_pathway(pathway_energies_temp)

    pickle.dump((
        rP.pickler,
        pathway_energies_temp,
        pathway_scores),
        open('/data/output_data/aziz_bcc_2'
			 '/pathway_energies_temp.p', 'wb'))
