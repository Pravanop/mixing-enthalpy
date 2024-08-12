"""This is the main script that will take an element list as input, the crystal system and other metadata to create
VASP runs for calculating energies of equimolar binary alloys.
Contributors : Joshua Cheng, Pravan Omprakash
"""
import os
from tqdm import tqdm

from mp_api.client import MPRester
from calculateEnthalpy.callMpAPI import eleList_to_POSCAR
from calculateEnthalpy.callMpAPI import getAPIKey
from getBinaryPairs.get_binary_pairs import get_binary_pairs
from obtainSQS.sqs_pymat import obtainSQS
from prepareVASPRuns.main import create_vasprun
from prepareVASPRuns.file_utils import load_yaml_to_dict

import warnings

warnings.filterwarnings("ignore")  # unsafe to use while testing.

main_input = load_yaml_to_dict("input.yaml")

benchmark = {
		'1' : True ,
		'2' : True ,
		'3' : True ,
		'4' : True ,
		'5' : True
		}
main_output_folder = f"{main_input['abs_path']}/Outputs_{main_input['mp-api']['lattice']}/"
if not os.path.exists(main_output_folder) :
	os.mkdir(main_output_folder)
"""
Step 1 - Get the lattice structures of all the elements in the input.
"""
out_file_path = f"{main_output_folder}/ele_POSCARs/"
if os.path.exists(out_file_path) :
	benchmark['1'] = False

if benchmark['1'] :
	mpr = MPRester(api_key = getAPIKey(f"{main_input['abs_path']}/callMpAPI/api_key.txt") , mute_progress_bars = True)
	ele_dict = eleList_to_POSCAR(
			mpr = mpr ,
			element_list = main_input['element_list'] ,
			lattice = main_input['mp-api']['lattice'] ,
			extract_style = "sym" ,
			out_file_path = out_file_path
			)
# raise ValueError("Stopping here")
"""
Step 2 - Get all binary pairs for the element list
"""

binary_pairs = get_binary_pairs(ele_list = main_input['element_list'])
print("Binary Pairs: ", binary_pairs)

"""
Step 3 - Create SQS for each pair
"""

sqs_folder_path = f"{main_output_folder}SQS_binaries"
if os.path.exists(sqs_folder_path) :
	benchmark['3'] = False

if benchmark['3'] :
	if not os.path.exists(sqs_folder_path) :
		os.mkdir(sqs_folder_path)
	
	for idx , ele_pair in enumerate(tqdm(binary_pairs , desc = "Creating SQS random alloy structures")) :
		
		ele_list = ele_pair.split('-')
		first_ele = ele_list[0]
		second_ele = ele_list[1]
		out_file_path2 = f"{sqs_folder_path}/{first_ele}{second_ele}"
		# if not os.path.exists(out_file_path2) :
		# 	os.mkdir(out_file_path2)
		inp_file_path = f"{out_file_path}/{first_ele}_{main_input['mp-api']['lattice']}.vasp"
		try:
			obtainSQS(
					supercell =  main_input['sqs']['supercell'],
					inp_file_path = inp_file_path ,
					output_path = out_file_path2 ,
					dopant_percentage = main_input['sqs']['doping_percent'] ,
					dopant = second_ele ,
					doping_site = first_ele ,
					)
		except:
			benchmark['3'] = False
			raise Exception("Something went wrong")
			

# """
# Step 4 Copying created POSCAR files to a new directory (I have a penchant for creating backups)
# """
#
# sqs_results_path = f"{main_output_folder}/SQS_results"
# if os.path.exists(sqs_results_path) :
# 	benchmark['4'] = False
#
# if benchmark['4'] :
# 	lsqs = os.listdir(sqs_folder_path)
#
# 	if not os.path.exists(sqs_results_path) :
# 		os.mkdir(sqs_results_path)
# 	for idx , dir in enumerate(tqdm(lsqs , desc = "Copying SQS POSCARs")) :
#
# 		shutil.copy(
# 				src = f"{sqs_folder_path}/{dir}.vasp" ,
# 				dst = f"{sqs_results_path}/{dir}.vasp"
# 				)

"""
Step 5 The final step of this half. Creating Vaspjob files.
"""

vasp_runs_path = f"{main_output_folder}/vasp_runs"
if os.path.exists(vasp_runs_path) :
	benchmark['5'] = False

if benchmark['5'] and main_input['run_vasp'] :
	lsqs_poscar = os.listdir(sqs_folder_path)
	
	if not os.path.exists(vasp_runs_path) :
		os.mkdir(vasp_runs_path)
	
	for idx , poscar in enumerate(tqdm(lsqs_poscar , desc = "Creating VASP runs")) :
		
		create_vasprun(
				inp_file_path = f"{sqs_folder_path}/{poscar}" ,
				output_folder_path = f"{vasp_runs_path}" ,
				metadata_path = f"{main_input['abs_path']}/vasp_metadata/metadata.yaml"
				)
