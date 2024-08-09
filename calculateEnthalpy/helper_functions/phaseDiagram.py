import json

import numpy as np
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core import Composition
from tqdm import tqdm

from calculateEnthalpy.helper_functions.data_utils import load_json
from calculateEnthalpy.helper_functions.grid_code import create_multinary
from calculateEnthalpy.helper_functions.thermo_math import thermo_maths
import warnings
warnings.filterwarnings("ignore")

class phaseDiagram:

    def __init__(self,
                 source: str = "pravan",
                 lattice: str = "bcc",
                 version_no: int = 1,
                 im_flag: bool = False,
                 ):

        if not im_flag:
            self.file_path = f"/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/data/output_data/{source}_{lattice}_{version_no}/multinaries_wo_im.json"
        else:
            self.file_path = f"/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/data/output_data/{source}_{lattice}_{version_no}/multinaries_im.json"

        with open(self.file_path, "r") as f:
            self.data = json.load(f)

        self.binary_dict = load_json(folder_path=f"/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/data/input_data/", lattice=lattice, source=source)
        self.tm = thermo_maths(self.binary_dict)
        self.temp_grid = np.arange(0, 3200, 200)

    def make_PD_composition(self, temperature: float, composition: str):
        """
        This should be a user-facing function in the future.
        :param composition:
        :param temperature:
        :return:
        """
        # comp is "El1-El2-EL3"
        global pd_entries_list
        self.ele_list = composition.split('-')
        n_alloy = len(self.ele_list)
        all_combs = create_multinary(element_list=self.ele_list, no_comb=list(range(2, n_alloy + 1)))
        pd_entry_input = {}
        for key, n_nary in all_combs.items():
            for idx2, subset_comp in enumerate(n_nary):
                temp_subset = self.data[str(key)][subset_comp]
                subset_list = subset_comp.split('-')

                for mol_ratio, value in temp_subset.items():
                    if mol_ratio == "intermetallic":
                        continue
                    mol_ratio = mol_ratio.split('-')
                    mol_ratio = np.array(mol_ratio).astype(float)
                    mol_ratio = dict(zip(subset_list, mol_ratio))
                    mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}
                    name = Composition(Composition(mol_ratio).get_integer_formula_and_factor()[0])
                    config_entropy = value['config_entropy']
                    pd_entry_input[name] = self.tm.calc_gibbs_energy(
                        value['mix_enthalpy'],
                        config_entropy,
                        temperature
                    ) * name.num_atoms

                # intermetallics
                if "intermetallic" in temp_subset:
                    for idx3, intermetallic in enumerate(temp_subset['intermetallic']):
                        name = Composition(intermetallic['formula_pretty'])
                        pd_entry_input[name] = intermetallic['formation_energy_per_atom'] * name.num_atoms

                # elements
                for ele in self.ele_list:
                    name = Composition(ele)
                    pd_entry_input[name] = 0 * name.num_atoms

            pd_entries_list = []
            for key, value in pd_entry_input.items():
                pd_entries_list.append(PDEntry(composition=key, energy=value))

        phase_diagram = PhaseDiagram(pd_entries_list)
        return phase_diagram

    def make_PD_comp_temp(self,
                          composition: str) -> dict[float:PhaseDiagram]:
        PD_temp_comp_dict = {}
        for idx, temp in enumerate(tqdm(self.temp_grid, desc="Running Temperature")):
            PD_temp_comp_dict[temp] = self.make_PD_composition(temperature=temp,
                                                               composition=composition)

        return PD_temp_comp_dict

    def check_stability(self, mol_ratio, temp, conv_hull):
        mol_ratio = dict(zip(self.ele_list, mol_ratio))
        mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}
        comp = '-'.join(list(mol_ratio.keys()))
        mix_enthalpy = self.tm.calc_multinary_mixEnthalpy(alloy_comp=comp, mol_ratio=mol_ratio)
        entropy = self.tm.calc_configEntropy(mol_ratio)
        gibbs = self.tm.calc_gibbs_energy(enthalpy=mix_enthalpy,
                                          entropy=entropy,
                                          temperature=temp)
        pdEntry = self._make_PD_entry(energy=gibbs, mol_ratio=mol_ratio)

        try:
            return conv_hull.get_decomp_and_e_above_hull(entry=pdEntry)

        except:
            return None

    def find_misc_temperature(self, mol_ratio, composition):
        for idx, temperature in enumerate(self.temp_grid):
            conv_hull = self.make_PD_composition(temperature=temperature,
                                                 composition=composition)
            is_stable = self.check_stability(mol_ratio=mol_ratio, temp=temperature, conv_hull=conv_hull)
            if is_stable is not None:
                if np.isclose(is_stable[1], 0.0, atol=1e-3):
                    return temperature

    def _make_PD_entry(self,
                       mol_ratio: dict,
                       energy: float):

        name = Composition(Composition(mol_ratio).get_integer_formula_and_factor()[0])
        return PDEntry(composition=name, energy=energy*name.num_atoms)
