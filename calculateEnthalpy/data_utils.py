"""
Helper functions.
"""
import json
import pandas as pd

def load_csv(
		folder_path: str = "./data/data_input" ,
		lattice: str = "bcc" ,
		source: str = "bokas"
		) -> pd.DataFrame :
	"""
	Loads a csv file from previous sources. Deprecated; json files are used instead.
	:param folder_path: the folder that holds all the input files from previous papers.
	:param lattice: values are 'bcc', 'fcc', 'hcp', 'sc'
	:param source: values for now are 'bokas', 'zhang', 'joshua'
	:return: a pandas dataframe
	"""
	file_path = f"{folder_path}/{lattice}_{source}.csv"
	df = pd.read_csv(file_path)
	df.set_index('ele' , inplace = True)
	return df

def load_json(
		folder_path: str = "./data/input_data" ,
		lattice: str = "bcc" ,
		source: str = "bokas"
		) -> dict :
	"""
	Loads a json file for mixing enthalpy values. MAIN METHOD TO USE!
	:param folder_path: the folder that holds all the input files from previous papers. TODO can make it more organized
	:param lattice: values are 'bcc', 'fcc', 'hcp', 'sc'
	:param source: values for now are 'bokas', 'zhang', 'joshua'
	:return: a dictionary
	"""
	with open(f'{folder_path}/{lattice}_{source}.json' , 'r') as f :
		return json.load(f)

def extract_ele_list(
		folder_path: str = "./data/input_data" ,
		lattice: str = "bcc" ,
		source: str = 'bokas'
		) -> list[str] :
	"""
	Removes the element list for a particular source. Needs to be stored with diligence in previous steps!
	:param source: values for now are 'bokas', 'zhang', 'joshua'
	:param lattice: values are 'bcc', 'fcc', 'hcp', 'sc'
	:param folder_path: the folder that holds all the input files from previous papers.
	:return: List of strings of element symbols
	"""
	with open(f"{folder_path}/element_list_{lattice}_{source}.txt" , 'r') as f :
		ele_list = f.read()
	
	return ele_list.split(',')

def dict_to_csv(
		results_dict: dict ,
		source: str ,
		lattice: str ,
		n_alloy: int ,
		output_folder: str = './data/data_output/'
		) -> None :
	"""

	:param results_dict: dump dictionary used for backup
	:param source: values for now are 'bokas', 'zhang', 'joshua'
	:param lattice: values are 'bcc', 'fcc', 'hcp', 'sc'
	:param output_folder:
	:param n_alloy:
	:return:
	"""
	
	useful_names = {
			2 : 'binary' ,
			3 : 'tertiary' ,
			4 : 'quaternary' ,
			5 : 'quinary' ,
			6 : 'senary' ,
			7 : 'septenary' ,
			8 : 'octanery'
			}
	
	temp_df = pd.DataFrame.from_dict(data = results_dict).transpose()
	path = f"{output_folder}/{useful_names[n_alloy]}_{lattice}_{source}_equiOut.csv"
	temp_df.to_csv(path)

def extract_binaryEnthalpy(
		binary_dict: dict ,
		ele_pair: str
		) -> float :
	"""
	extracts binary_enthalpy for element pair
	:param binary_dict: The required dict from the json file
	:param ele_pair: a str of the form 'E1-E2'
	:return: binary enthalpy value in eV/atom
	"""
	return binary_dict[ele_pair]