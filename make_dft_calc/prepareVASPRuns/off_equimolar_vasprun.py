from make_dft_calc.prepareVASPRuns.main import create_vasprun
import os
from tqdm import tqdm

system = "TiHf"
lattice = "BCC"
sqs_results_path = f"Outputs_off_equimolar_{lattice}/{system}"
vasp_runs_path = f"Outputs_off_equimolar_{lattice}/{system}/{system}_{lattice}_vasp"

lsqs_poscar = os.listdir(sqs_results_path)

if not os.path.exists(vasp_runs_path):
	os.mkdir(vasp_runs_path)

for idx, poscar in enumerate(tqdm(lsqs_poscar, desc="Creating VASP runs")):
	create_vasprun(
		inp_file_path=f"{sqs_results_path}/{poscar}",
		output_folder_path=f"{vasp_runs_path}",
		metadata_path=f"/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/vasp_metadata/metadata.yaml"
	)
