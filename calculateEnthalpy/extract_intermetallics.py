"""
Calls MP-API to get all intermetallics for everything from previous step.

"""

import json
from data_utils import load_json
from callMpAPI.utils import getAPIKey
from mp_api.client import MPRester
import time

def extract_intermetallic(comp , fields , mpr) :
	docs = mpr.materials.summary.search(chemsys = comp , fields = fields)
	# print(docs)
	return docs

def intermetallic_sorter(n_nary_dict , docs) :
	
	for idx , doc in enumerate(docs) :
		
		temp_dict = {
				'formula_pretty'            : str(doc.composition.reduced_formula) ,
				'formation_energy_per_atom' : doc.formation_energy_per_atom ,
				'energy_above_hull'         : doc.energy_above_hull
				}
		comp = str(doc.chemsys)
		if "intermetallic" in n_nary_dict[comp] :
			n_nary_dict[comp]['intermetallic'].append(temp_dict)
		else :
			n_nary_dict[comp]['intermetallic'] = [temp_dict]
	
	return n_nary_dict

def find_n_nary_intermetallic(n_nary_dict , mpr , fields) :
	
	comps = list(n_nary_dict.keys())
	docs = extract_intermetallic(
		comp = comps ,
		mpr = mpr ,
		fields = fields
		)
	n_nary_dict = intermetallic_sorter(
		n_nary_dict = n_nary_dict ,
		docs = docs
		)
	
	return n_nary_dict

start = time.time()

def extract_intermetallics():
	"""
	
	I have put all the script code into this function. When a final script is made to put all this together, then we can
	work on making it compatible, but looks like only a few parameters need to be updated
	:return:
	"""
	
	mpr = MPRester(api_key = getAPIKey("/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/callMpAPI/api_key.txt"))
	fields = ['composition' , 'formation_energy_per_atom' , 'energy_above_hull' , 'chemsys']
	
	with open(
            "/calculateEnthalpy/data/output_data/old_data/dump_20240404-213329.json", 'r'
			) as f :
		dump_dict = json.load(f)
	
	ternary = dump_dict['3']
	quartenary = dump_dict['4']
	quinary = dump_dict['5']
	
	# load binary
	binary = load_json(source = "Cr_W", lattice = "bcc")  # taking default values for now. #BCC, bokas, input_data/
	
	binary_new = {}
	for key , value in binary.items() :
		binary_new[key] = {
				'mix_enthalpy'   : value ,
				'config_entropy' : 5.97e-5 ,
				}
	
	ternary_updated = find_n_nary_intermetallic(ternary , mpr , fields)
	binary_updated = find_n_nary_intermetallic(binary_new , mpr , fields)
	stop = time.time()
	print(stop - start)
	
	dump_dict['2'] = binary_updated
	dump_dict['3'] = ternary_updated
	dump_dict['4'] = quartenary
	dump_dict['5'] = quinary
	
	with open(f'data/output_data/dump_{time.strftime("%Y%m%d-%H%M%S")}_intermetallic.json' , 'w') as f :
		json.dump(dump_dict , f , ensure_ascii = False , indent = 4)


extract_intermetallics()