import itertools
import json
from typing import Union

import pandas as pd
from pandas import DataFrame
from calculateEnthalpy.helper_functions.grid_code import create_multinary
from calculateEnthalpy.helper_functions.phaseDiagram import phaseDiagram

# open intermetallics file


# We need to define standard input format. Best is to say ['Al', 'W', 'Ni']

im_path = "calculateEnthalpy/data/output_data/pravan_bcc_1/multinaries_wo_im.json"


class search_composition:
    def __init__(
            self,
            inp: list,
            path: str = im_path,
    ):
        inp.sort()
        self.ele_list = inp
        self.primary_key = str(len(self.ele_list))
        self.secondary_key = "-".join(self.ele_list)
        self.n_nary_map = {
            'binary': 2,
            'ternary': 3,
            'quaternary': 4,
            'quinary': 5
        }

        with open(path, 'r') as f:
            self.intermetallic_dict = json.load(f)

    def search_comp(self, ele_list):
        primary_key = str(len(ele_list))
        secondary_key = "-".join(ele_list)
        try:
            result = self.intermetallic_dict[primary_key][secondary_key]
            result.pop('intermetallic', None)
            return {'Composition': ''.join(secondary_key.split('-')), 'Enthalpy': result['mix_enthalpy']}
        except KeyError as e:
            return "One of your elements is not in the database, please check the list below."

    @property
    def search_equimolar(self) -> str:
        """

		"""
        if len(self.ele_list) > 5:
            return "Sorry! We currently do not have compositions greater than quinary!"
        answer = self.search_comp(self.ele_list)
        are = ''
        for key, value in answer.items():
            are += str(key) + ' : ' + str(value) + " , "
        return are

    @property
    def search_all_intermetallics(self) -> Union[pd.DataFrame, str]:

        if len(self.ele_list) > 5:
            return "Sorry! We currently do not have compositions greater than quinary!"
        try:
            result = self.intermetallic_dict[self.primary_key][self.secondary_key]['intermetallic']
            result_df = pd.DataFrame.from_dict(result, orient='columns')
            return result_df
        except KeyError as e:
            return "There are no intermetallics from materials project database"

    def search_combinations(self, combinations: list[str]) -> Union[str, DataFrame]:
        """

		:return:
		"""
        self.combinations = combinations
        while None in self.combinations:
            self.combinations.remove(None)

        if self.combinations:
            if self.combinations is str:
                no_combs = [self.combinations]
            else:
                no_combs = [self.n_nary_map[i] for i in self.combinations]

            if any(y > int(self.primary_key) for y in no_combs):
                greater = [i for i in no_combs if i > int(self.primary_key)]
                [no_combs.remove(i) for i in greater]
                if not no_combs:
                    return "You have asked for an illegal combination. Please try again."

            combs = create_multinary(self.ele_list, no_comb=no_combs)
            total_combs = list(itertools.chain.from_iterable(combs))
            ele_lists = [i.split('-') for i in total_combs]
            answer = []
            for ele_list in ele_lists:
                answer.append(self.search_comp(ele_list))

            df = pd.DataFrame(answer)
            if df.empty:
                return "Could not find the combination(s) you were looking for."
            else:
                return df
        else:
            return "Please pick a combination!"

    def get_phase_diagram(self, temperature, flag="offequi"):
        # print(self.secondary_key)
        if flag == "equi":
            answer = phaseDiagram(composition=self.secondary_key).make_PD_composition(temperature)
        if flag == "offequi":
            answer = phaseDiagram(composition=self.secondary_key).make_PD_composition(temperature)
        return answer


class search_any:
    def __init__(
            self,
            inp: Union[list, str],
            find: list[str],
            path: str = im_path,
            element_list_path: str = "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/data/input_data/pravan/element_list_bcc_pravan.txt",
    ):

        with open(path, 'r') as f:
            self.intermetallic_dict = json.load(f)
        with open(element_list_path, 'r') as f:
            self.element_total_list = f.read()

        if isinstance(inp, str):
            self.ele = inp
            inp = [inp]

        if isinstance(inp, list):
            inp.sort()
            self.ele = '-'.join(inp)

        self.n_nary_map = {
            'binaries': 2,
            'ternaries': 3,
            'quaternaries': 4,
            'quinaries': 5
        }

        self.find = find
        self.process_find()

    def process_find(self):
        self.find_processed = [self.n_nary_map[i] for i in self.find]

    @property
    def search_all_contains(self) -> Union[list, str]:
        """

		:return:
		"""
        if len(self.ele.split('-')) > 5:
            return "Sorry! We currently do not have compositions greater than quinary!"
        values = []
        for primary_key in self.find_processed:
            temp_dict = self.intermetallic_dict[str(primary_key)]

            for key, value in temp_dict.items():
                if self.ele in key:
                    values.append({'Composition': key, 'Enthalpy': value['mix_enthalpy']})

        return values

# search_inst = search_any(inp = ['Al', 'Ni', 'Fe'], find = ['binary', 'ternary', 'quaternary'])
