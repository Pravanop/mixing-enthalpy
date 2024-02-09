import warnings
warnings.filterwarnings("ignore" , category = DeprecationWarning)
from typing import Union
from itertools import combinations
import json
import numpy as np
import pandas as pd
from time import time
import datetime
from tqdm import tqdm



"""
I really liked your error handling, but I think we can make it very clean now.
"""

"""
The input file is an excel file. Best we keep it as a separate csv, that is untouched. Regardless, we have a 2D grid
of enthalpy values.
"""
"""
The user can do multiple things on this 2D grid.
They can want a specific entry. We need a function for that. I think that is done.
Secondly they can provide a group of elements and ask for all possible combinations (we ourselves need this many times).
Thirdly they can just ask for one element(s), and extract all binaries, ternaries etc.
We need three input use cases for this.
"""
"""
Apart from that, the main functionality is: taking a data point(s), then performing some maths on it (this is my
naive thinking). I need to sit down with you to understand the maths, and we need to clearly document it and store it away, never to be touched again.
The other part which seems to take a lot of your code up is the screening of good alloys, and creation of phase
diagrams. I am not yet able to understand how that works, so I will leave it for later, I am going to copy/modify
your code for the first part herein.
"""

"""
Helper functions.
"""

def load_csv(
		folder_path: str = "data" ,
		lattice: str = "bcc" ,
		source: str = "bokas"
		) :
	"""
	
	:param folder_path:
	:param lattice:
	:param source:
	:return:
	"""
	file_path = f"{folder_path}/{lattice}_{source}.csv"
	df = pd.read_csv(file_path)
	# df.set_index('ele')
	return df

def load_json(
		folder_path: str = "data" ,
		lattice: str = "bcc" ,
		source: str = "bokas"
		) :
	"""
	
	:param folder_path:
	:param lattice:
	:param source:
	:return:
	"""
	with open('data/input_data/bcc_bokas.json' , 'r') as f :
		return json.load(f)

def extract_ele_list(
		file_path: str = "./data/input_data" ,
		lattice: str = "bcc" ,
		source: str = 'bokas'
		) :
	"""
	
	:param df:
	:return:
	"""
	with open(f"{file_path}/element_list_{lattice}_{source}.txt" , 'r') as f :
		ele_list = f.read()
	
	return ele_list.split(',')

def dict_to_csv(
		results_dict: dict ,
		source: str ,
		lattice: str ,
		output_folder: str ,
		n_alloy: int
		) :
	
	useful_names = {
			2 : 'binary' ,
			3 : 'tertiary' ,
			4 : 'quartenary' ,
			5 : 'quinary' ,
			6 : 'senary' ,
			7 : 'septenary' ,
			8 : 'octanery'
			}
	
	temp_df = pd.DataFrame.from_dict(data = results_dict).transpose()
	path = f"{output_folder}/{useful_names[n_alloy]}_{lattice}_{source}_out.csv"
	temp_df.to_csv(path)

"""
The main logic
"""

def extract_binaryEnthalpy(
		binary_dict: dict ,
		ele_pair: str
		) :
	"""
	
	:param ele_pair:
	:return:
	"""
	return binary_dict[ele_pair]

def create_multinary(
		element_list: list[str] ,
		no_comb: Union[list[int] , str]
		) :
	"""
	I am just putting together all the functions you had written into one function, because they all do the same thing.
	:param element_list: The list of element symbols to create combinations for
	:param no_comb: A list of number to create combinations for. For ex. [2,3] means create binary, ternary. Can take value "all" for combinations from 2-7.
	:return:
	"""
	if no_comb == "all" :
		no_comb = list(range(2 , 8))
	
	result = []
	for idx , comb in enumerate(no_comb) :
		temp_joint = list(combinations(element_list , comb))
		s = '-'  # To join the element
		final_joint = sorted(
				{s.join(sorted(pair)) for pair in temp_joint}
				)
		if final_joint != [] :  # incase the multinary asked is greater than list length
			result.append(final_joint)
	
	return result

"""Maths"""

def calc_pairwiseInteractionParameter(mix_enthalpy: float , mol_i: float , mol_j: float) :
	"""
	Function for the eq 1 in Zhang et al. A Fast and Robust Method for Predicting the Phase Stability of Refractory
	Complex Concentrated Alloys Using Pairwise Mixing Enthalpy. http://dx.doi.org/10.2139/ssrn.4081906
	delta_H_ij = omega_ij*x_i*x_j
	
	:param mix_enthalpy:
	:param mol_i:
	:param mol_j:
	:return:
	"""
	return mix_enthalpy / (mol_i * mol_j)

def calc_multinary_mixEnthalpy(
		alloy_comp: str ,
		binary_dict: dict ,
		mol_ratio: float
		) :
	"""
	Function for eq 2 in Zhang et al. A Fast and Robust Method for Predicting the Phase Stability of Refractory
	Complex Concentrated Alloys Using Pairwise Mixing Enthalpy. http://dx.doi.org/10.2139/ssrn.4081906
	delta_Hmix = sum_i=1,i!=j_n(omega_ij*x_i*x_j)
	:param alloy_comp: str of the form 'ele1-ele2-ele3'
	:param df:
	:param mol_ratio:
	:return:
	"""
	ele_list = alloy_comp.split('-')
	binaries = create_multinary(element_list = ele_list , no_comb = [2])[0]
	mix_enthalpy = 0
	for idx , binary in enumerate(binaries) :
		delH_ij = extract_binaryEnthalpy(
				binary_dict = binary_dict ,
				ele_pair = binary
				)
		omega_ij = calc_pairwiseInteractionParameter(
				mix_enthalpy = delH_ij ,
				mol_i = mol_ratio ,
				mol_j = 1 - mol_ratio
				)
		mix_enthalpy += omega_ij * mol_ratio * (1 - mol_ratio)
	
	return mix_enthalpy

def calc_configEntropy(mol_ratio: list[float]) :
	"""
	A simple function to calculate boltzmann configurational entropy.
	delta_S = -k_b*sum_i=1_n(x_i*ln(x_i))
	:param alloy_comp:
	:param mol_ratio:
	:return:
	"""
	assert sum(mol_ratio) == 1
	boltzmann_constant = 1e-1  # TODO change this to actual number
	return np.round(-boltzmann_constant * sum([x * np.log(x) for x in mol_ratio]), 4)

"""Main function for our codes"""

def calc_mixEnthalpy_dataset(
		folder_path: str = "./input_data/data" ,
		lattice: str = "bcc" ,
		source: str = "bokas" ,
		mol_ratio: float = 0.5 ,
		max_alloy_n: int = 5
		) :
	"""
	This function ties everything up together. I hope this is an efficient way. Subject to change.
	:param folder_path:
	:param lattice:
	:param source:
	:param mol_ratio:
	:param max_alloy_n: maximum size alloy. 5 means quinary, 3 means tertiary etc.
	:return:
	"""
	start = time()
	print("==========")
	print("Status Update")
	print("==========")
	print(f"Loading binary enthalpies of {lattice} phase from {source} et al.")
	# let us get the binary file
	binary_dict = load_json(folder_path = folder_path , lattice = lattice , source = source)
	# access the element list
	ele_list = extract_ele_list()
	print("==========")
	print("Status Update")
	print("==========")
	print(f"Creating 3-nary to {max_alloy_n}-nary equimolar alloy compositions from {len(ele_list)} elements")
	
	# let us create all the pairs possible upto 5 for now
	multinaries = create_multinary(element_list = ele_list , no_comb = list(range(3 , max_alloy_n + 1)))
	# now we have the multinary_pairs, we can start to calculate the enthalpies and store each one neatly.
	# let us begin with ternary and work our way up
	results_dict = {}
	print("==========")
	print("Status Update")
	print("==========")
	
	
	for idx , multinary in enumerate(multinaries) :
		# take each alloy in this now
		temp_zero = multinary[0]
		n_alloy = len(temp_zero.split('-'))  # just to calculate the number of elements
		multinary_dict = {}
		for alloy_idx , alloy in enumerate(tqdm(multinary, desc = f"Accessing {n_alloy}-nary")) :
			
			mix_enthalpy = calc_multinary_mixEnthalpy(
					alloy_comp = alloy ,
					binary_dict = binary_dict ,
					mol_ratio = 0.5
					)
			equimol_ratio = [1 / n_alloy] * n_alloy
			entropy = calc_configEntropy(mol_ratio = equimol_ratio)
			multinary_dict[alloy] = {
					'mix_enthalpy'   : mix_enthalpy ,
					'config_entropy' : entropy
					}
		
		dict_to_csv(
				results_dict = multinary_dict ,
				source = source ,
				lattice = lattice ,
				output_folder = './data/output_data/' , #ideally should not take her
				n_alloy = n_alloy
				)
		
		results_dict[n_alloy] = multinary_dict  # so this should create a dict(dict(dict)))
		
	# Backup, backup, backup the backup.
	with open(f'data/output_data/dump.json' , 'w') as f :  # TODO add datetime
		json.dump(results_dict , f , ensure_ascii = False , indent = 4)
	stop = time()
	print("==========")
	print(f"Code has run upto {max_alloy_n}-nary alloy for {len(ele_list)} elements in {np.round(stop - start, 3)} "
	      f"s.")
	print("==========")
	return results_dict

calc_mixEnthalpy_dataset()
