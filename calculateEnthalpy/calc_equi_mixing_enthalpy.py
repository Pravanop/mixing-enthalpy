"""
Calculating mixing enthalpy of any n-nary alloy using pairwise equimolar binary mixing enthalpies.
Contributors: Joshua Cheng, Pravan Omprakash (M-CUBE)

"""
import warnings

warnings.filterwarnings("ignore" , category = DeprecationWarning)
import json
import numpy as np
import pandas as pd
import time

from typing import Union , Any
from itertools import combinations, permutations
from tqdm import tqdm

from data_utils import *
from thermo_math import *
from create_alloy_comp import *

"""
I really liked your error handling, but I think we can make it very clean now, and then you can add those parts later.
"""


def calc_mixEnthalpy_dataset(
		folder_path: str = "./data/input_data" ,
		lattice: str = "bcc" ,
		source: str = "bokas" ,
		max_alloy_n: int = 5
		) -> dict[int , dict[Any , dict[str , float]]] :
	"""
	This function creates n-nary equimolar compositions
	The steps I have followed:
	1) Extract binary enthalpy values from a specific source.
	2) Make combinations of element list from a specific source up to a specified n-nary alloy.
	3) Calculate mixing enthalpy and entropy for n-nary alloys as in the paper by Zhang et al. (found in docs).
	4) Created separate csv files for each n-nary alloy and also backup files. Everything is found in data folder.
	
	:param folder_path: the folder that holds all the input files from previous papers.
	:param lattice: values are 'bcc', 'fcc', 'hcp', 'sc'
	:param source: values for now are 'bokas', 'zhang', 'joshua'
	:param mol_ratio: binary mol_ratio
	:param max_alloy_n: maximum size alloy. 5 means quinary, 3 means tertiary etc.
	:return: a results dict that has been dumped
	"""
	
	start = time.time()
	
	# Step 1
	print("==========")
	print("Status Update")
	print("==========")
	print(f"Loading binary enthalpies of {lattice} phase from {source} et al.")
	
	binary_dict = load_json(folder_path = folder_path , lattice = lattice , source = source)
	# access the element list
	ele_list = extract_ele_list(folder_path=folder_path, lattice=lattice, source=source) # fixed bug
	
	# Step 2
	print("==========")
	print("Status Update")
	print("==========")
	print(f"Creating 3-nary to {max_alloy_n}-nary equimolar alloy compositions from {len(ele_list)} elements")
	
	multinaries = create_multinary(element_list = ele_list , no_comb = list(range(3 , max_alloy_n + 1)))
	
	# Step 3
	print("==========")
	print("Status Update")
	print("==========")
	non_equi_binary_dict = load_json(folder_path = folder_path , lattice = lattice , source = source + "_offEqui")
	results_dict = {}
	for idx , multinary in enumerate(multinaries) :
		# take each alloy in this now
		temp_zero = multinary[0]
		n_alloy = len(temp_zero.split('-'))  # just to calculate the number of elements
		multinary_dict = {}
		for alloy_idx , alloy in enumerate(tqdm(multinary , desc = f"Accessing {n_alloy}-nary")) :
			ele_list = alloy.split('-')
			equimol_ratio = [1 / n_alloy] * n_alloy
			mol_ratio = dict(zip(ele_list, equimol_ratio))
			try:
				equi_mix_enthalpy = calc_multinary_mixEnthalpy(
						alloy_comp = alloy ,
						binary_dict = binary_dict ,
						mol_ratio = mol_ratio
						)
			except KeyError as e:
				continue
			entropy = calc_configEntropy(mol_ratio = mol_ratio)

			x_mol = np.linspace(0, 1, 10)[1:]
			comb_mol = np.array(list(permutations(x_mol , n_alloy)))
			off_equi = {}
			for mol_idx, mol in enumerate(comb_mol):
				
				if np.sum(mol) == 1:
					mol_ratio = dict(zip(ele_list , mol))
					
					try :
						non_equi_mix_enthalpy = calc_multinary_mixEnthalpy_offEqui(
								alloy_comp = alloy ,
								binary_dict = non_equi_binary_dict ,
								mol_ratio = mol_ratio
								)
						
					except KeyError as e :
						break
					off_equi[','.join(np.array(mol).astype(str))] = non_equi_mix_enthalpy
					
				else:
					continue
				
			multinary_dict[alloy] = {
					'off_equi' : off_equi,
					'mix_enthalpy'   : equi_mix_enthalpy ,
					'config_entropy' : entropy
					}
		
		dict_to_csv(
				results_dict = multinary_dict ,
				source = source ,
				lattice = lattice ,
				output_folder = './data/output_data/' ,  # ideally should not take here
				n_alloy = n_alloy
				)
		
		results_dict[n_alloy] = multinary_dict
	
	# Backup, backup, backup the backup.
	with open(f'data/output_data/dump_{time.strftime("%Y%m%d-%H%M%S")}.json' , 'w') as f :
		json.dump(results_dict , f , ensure_ascii = False , indent = 4)
	
	print("==========")
	print(f"The output file has been dumped with timestamp at data/output_data/")
	
	stop = time.time()
	
	print("==========")
	print(
			f"Code has run up to {max_alloy_n}-nary alloy for {len(ele_list)} elements in "
			f"{np.round(stop - start , 3)} "
			f"s."
			)
	print("==========")
	return results_dict


calc_mixEnthalpy_dataset(lattice = "bcc", source = "Cr_W")
