"""This is the main script that will take an element list as input, the crystal system and other metadata to create
VASP runs for calculating energies of equimolar binary alloys.
Contributors : Joshua Cheng, Pravan Omprakash
"""
import os
from tqdm import tqdm
import shutil

from mp_api.client import MPRester
from callMpAPI.extract import eleList_to_POSCAR
from callMpAPI.utils import getAPIKey
from getBinaryPairs.get_binary_pairs import get_binary_pairs
from obtainSQS.main import write_rndstr , write_sqscellout
from obtainSQS.run_atat import runsqs
from obtainSQS.sqs2poscar import sqs2POSCAR
from prepareVASPRuns.main import create_vasprun
from prepareVASPRuns.file_utils import load_yaml_to_dict


main_input = load_yaml_to_dict("input.yaml")

benchmark = {
		'1': True,
		'2': True,
		'3': True,
		'4': True,
		'5': True
		}
main_output_folder = f"{main_input['abs_path']}/Outputs/"
if not os.path.exists(main_output_folder) :
	os.mkdir(main_output_folder)
"""
Step 1 - Get the lattice structures of all the elements in the input.
"""
out_file_path = f"{main_input['abs_path']}/Outputs/ele_POSCARs/"
if os.path.exists(out_file_path):
	benchmark['1'] = False
	
if benchmark['1']:
	mpr = MPRester(api_key = getAPIKey(f"{main_input['abs_path']}/callMpAPI/api_key.txt") , mute_progress_bars = True)
	ele_dict = eleList_to_POSCAR(
			mpr = mpr ,
			element_list = main_input['element_list'] ,
			lattice = main_input['mp-api']['lattice'] ,
			extract_style = "sym" ,
			out_file_path = out_file_path
			)

"""
Step 2 - Get all binary pairs for the element list
"""

binary_pairs = get_binary_pairs(ele_list = main_input['element_list'])

"""
Step 3 - Create SQS for each pair
"""

sqs_folder_path = f"{main_input['abs_path']}/Outputs/SQS_binaries"
if os.path.exists(sqs_folder_path):
	benchmark['3'] = False
	
if benchmark['3']:
	if not os.path.exists(sqs_folder_path) :
		os.mkdir(sqs_folder_path)
	
	for idx , ele_pair in enumerate(tqdm(binary_pairs , desc = "Creating SQS random alloy structures")) :
		
		ele_list = ele_pair.split('-')
		first_ele = ele_list[0]
		second_ele = ele_list[1]
		out_file_path = f"{sqs_folder_path}/{first_ele}{second_ele}/"
		if not os.path.exists(out_file_path) :
			os.mkdir(out_file_path)
		
		structure = ele_dict[first_ele]
		write_rndstr(
				structure = structure ,
				output_path = out_file_path ,
				doping_percent = main_input['sqs']['doping_percent'] ,
				doping_site = first_ele ,
				doping_ele = second_ele,
				coord_type = main_input['sqs']['coord_type']
				)
		
		write_sqscellout(
				supercell = main_input['sqs']['supercell'] ,
				output_path = out_file_path
				)
		
		sqs_status = runsqs(
				corrdump_path = main_input['sqs']['corrdump_path'],
				mcsqs_path = main_input['sqs']['mcsqs_path'] ,
				path_dir = out_file_path ,
				pair_interaction = main_input['sqs']['pair_interaction'] ,
				triplet_interaction = main_input['sqs']['triplet_interaction']
				)  # TODO make nn_dist functions
		
		if sqs_status :
			sqs2POSCAR(
					input_file_path = f"{out_file_path}/bestsqs.out" ,
					output_file_path = f"{out_file_path}/{first_ele}{second_ele}"
					)
		
		else :
			print("Error")  # TODO raise exceptions

"""
Step 4 Copying created POSCAR files to a new directory (I have a penchant for creating backups)
"""

sqs_results_path = f"{main_input['abs_path']}/Outputs/SQS_results"
if os.path.exists(sqs_results_path):
	benchmark['4'] = False

if benchmark['4']:
	lsqs = os.listdir(sqs_folder_path)
	
	if not os.path.exists(sqs_results_path) :
		os.mkdir(sqs_results_path)
	for idx , dir in enumerate(tqdm(lsqs , desc = "Copying SQS POSCARs")) :
		
		shutil.copy(
			src = f"{sqs_folder_path}/{dir}/{dir}.vasp" ,
			dst = f"{sqs_results_path}/{dir}.vasp"
			)

"""
Step 5 The final step of this half. Creating Vaspjob files.
"""

vasp_runs_path = f"{main_input['abs_path']}/Outputs/vasp_runs"
if os.path.exists(vasp_runs_path):
	benchmark['5'] = False

if benchmark['5'] and main_input['run_vasp']:
	lsqs_poscar = os.listdir(sqs_results_path)
	
	if not os.path.exists(vasp_runs_path) :
		os.mkdir(vasp_runs_path)
	
	for idx , poscar in enumerate(tqdm(lsqs_poscar , desc = "Creating VASP runs")) :
		
		create_vasprun(
				inp_file_path = f"{sqs_results_path}/{poscar}" ,
				output_folder_path = f"{vasp_runs_path}" ,
				metadata_path = f"{main_input['abs_path']}/vasp_metadata/metadata.yaml"
				)
