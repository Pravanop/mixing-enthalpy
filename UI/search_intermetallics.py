import itertools
import json
from typing import Union

import pandas as pd
from pandas import DataFrame
from calculateEnthalpy.phase_diagram import make_PD_per_comp
from calculateEnthalpy.create_alloy_comp import create_multinary

# open intermetallics file


# We need to define standard input format. Best is to say ['Al', 'W', 'Ni']

class search() :
	def __init__(
			self , inp: list ,
			path: str = "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/data/output_data"
			            "/dump_20240217"
			            "-165153_intermetallic.json" ,
			**kwargs
			) :
		inp.sort()
		self.ele_list = inp
		self.primary_key = str(len(self.ele_list))
		self.secondary_key = "-".join(self.ele_list)
		self.n_nary_map = {
				'binary'     : 2 ,
				'ternary'    : 3 ,
				'quartenary' : 4 ,
				'quinary'    : 5
				}
		
		if kwargs['combinations'] :
			self.combinations = kwargs['combinations']
		
		with open(path , 'r') as f :
			self.intermetallic_dict = json.load(f)
	
	def search_comp(self , ele_list) :
		primary_key = str(len(ele_list))
		secondary_key = "-".join(ele_list)
		try :
			result = self.intermetallic_dict[primary_key][secondary_key]
			result.pop('intermetallic' , None)
			return {'Composition' : ''.join(secondary_key.split('-')) , 'Enthalpy' : result['mix_enthalpy']}
		except KeyError as e :
			return "One of your elements is not in the database, please check the list below."
	
	@property
	def search_equimolar(self) -> str :
		"""

		"""
		answer = self.search_comp(self.ele_list)
		are = ''
		for key , value in answer.items() :
			are += str(key) + ': ' + str(value) + ","
		return are
	
	@property
	def search_all_intermetallics(self) -> Union[pd.DataFrame , str] :
		try :
			result = self.intermetallic_dict[self.primary_key][self.secondary_key]['intermetallic']
			result_df = pd.DataFrame.from_dict(result , orient = 'columns')
			return result_df
		except KeyError as e :
			return "There are no intermetallics from materials project database"
	
	def search_combinations(self) -> Union[str , DataFrame] :
		"""

		:return:
		"""
		if self.combinations is str :
			no_combs = [self.combinations]
		else:
			no_combs = [self.n_nary_map[i] for i in self.combinations]
		
		if any(y > int(self.primary_key) for y in no_combs) :
			greater = [i for i in no_combs if i > int(self.primary_key)]
			[no_combs.remove(i) for i in greater]
			print("You have asked for an illegal combination. Please try again.")
		
		combs = create_multinary(self.ele_list , no_comb = no_combs)
		total_combs = list(itertools.chain.from_iterable(combs))
		ele_lists = [i.split('-') for i in total_combs]
		answer = []
		for ele_list in ele_lists :
			answer.append(self.search_comp(ele_list))
		
		df = pd.DataFrame(answer)
		if df.empty :
			return "Could not find the combination(s) you were looking for."
		else :
			return df
	
	def get_phase_diagram(self, temperature):
		
		answer = make_PD_per_comp(self.secondary_key, self.intermetallic_dict, temperature)
		return answer

search_inst = search(inp = ['Cr' , 'Mo' , 'Ta' , 'Pd'] , combinations = ['binary' , 'ternary'])

print(search_inst.get_phase_diagram(temperature = 1000))
