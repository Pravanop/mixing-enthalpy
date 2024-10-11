from make_dft_calc.prepareVASPRuns.main import create_vasprun
import os
from tqdm import tqdm

lattice = "BCC"
sqs_results_path = f"../Outputs_{lattice}"
# sqs_results_path = f"../../elements/{lattice}"

lfolder = os.listdir(sqs_results_path)
vasp_runs_path = f"../Outputs_{lattice}_runs/"
for i in tqdm(lfolder, desc = "making VASP runs"):
	temp_lfolder = os.listdir(f"{sqs_results_path}/{i}")
	output_path = vasp_runs_path + f"{i}"
	if not os.path.exists(output_path):
		os.mkdir(output_path)
	for poscar in temp_lfolder:
		create_vasprun(
			inp_file_path=f"{sqs_results_path}/{i}/{poscar}",
			output_folder_path=f"{output_path}",
			metadata_path=f"/Users/mcube/Desktop/Projects/make_dft_calc/prepareVASPRuns/vasp_metadata/metadata.yaml"
		)

# for i in tqdm(lfolder, desc = "making VASP runs"):
# 	vasp_runs_path = f"../../elements/{lattice}_runs"
# 	create_vasprun(
# 			inp_file_path=f"{sqs_results_path}/{i}",
# 			output_folder_path=f"{vasp_runs_path}",
# 			metadata_path=f"/Users/mcube/Desktop/Projects/make_dft_calc/prepareVASPRuns/vasp_metadata/metadata.yaml"
# 		)