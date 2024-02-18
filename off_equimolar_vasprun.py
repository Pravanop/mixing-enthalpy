from prepareVASPRuns.main import create_vasprun
import os
from tqdm import tqdm

sqs_results_path = "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/Outputs_off_equimolar/TaW"
vasp_runs_path = "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/Outputs_off_equimolar/TaW_vasp"

lsqs_poscar = os.listdir(sqs_results_path)

if not os.path.exists(vasp_runs_path) :
	os.mkdir(vasp_runs_path)

for idx , poscar in enumerate(tqdm(lsqs_poscar , desc = "Creating VASP runs")) :
	
	create_vasprun(
			inp_file_path = f"{sqs_results_path}/{poscar}" ,
			output_folder_path = f"{vasp_runs_path}" ,
			metadata_path = f"/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/vasp_metadata/metadata.yaml"
			)