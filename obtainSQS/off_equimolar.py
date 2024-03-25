from pymatgen.core.structure import Structure
from obtainSQS.main import write_rndstr , write_sqscellout
from obtainSQS.sqs2poscar import sqs2POSCAR
from obtainSQS.run_atat import runsqs
import os
from tqdm import tqdm
import shutil

def obtainSQS_off_equimolar(
		file_path: str ,
		system: str ,
		elements: list ,
		out_folder_path: str
		) :
	"""

	:param file_path:
	:param system:
	:param elements:
	:param out_folder_path:
	"""
	structure = Structure.from_file(f"{file_path}")
	
	for doping_percent in [12.5 , 25 , 37.5 , 50 , 62.5 , 75 , 87.5] :
		project_name = f"{system}_{doping_percent}"
		
		out_file_path = f"{out_folder_path}{project_name}"
		if not os.path.exists(out_file_path) :
			os.mkdir(out_file_path)
		
		write_rndstr(
				structure = structure ,
				output_path = out_file_path ,
				doping_percent = doping_percent ,
				doping_site = elements[0] ,
				doping_ele = elements[1] ,
				coord_type = 'frac'
				)
		
		write_sqscellout(
				supercell = [3 , 2 , 2] ,
				output_path = out_file_path
				)
		
		sqs_status = runsqs(
				corrdump_path = "/Users/pravanomprakash/Documents/Projects/atat/corrdump" ,
				mcsqs_path = "/Users/pravanomprakash/Documents/Projects/atat/mcsqs" ,
				path_dir = out_file_path ,
				pair_interaction = 2.4 ,
				triplet_interaction = 3.0
				)
		
		sqs2POSCAR(
				input_file_path = f"{out_file_path}/bestsqs.out" ,
				output_file_path = f"{out_file_path}/{project_name}"
				)


start = "Ti"
end = "Hf"
lattice = "BCC"
file_path = f"/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/Outputs_BCC/ele_POSCARs/{start}_{lattice}.vasp"
out_file_path = "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/Outputs_off_equimolar_BCC/"
system = f"{start}{end}"
os.mkdir(f"{out_file_path}{system}")
element_list = [start , end]
obtainSQS_off_equimolar(
		file_path = file_path ,
		out_folder_path = f"{out_file_path}{system}/" ,
		system = system ,
		elements = element_list
		)

lsqs = os.listdir(f"{out_file_path}{system}")
if not os.path.exists(f"{out_file_path}{system}_poscar/") :
	os.mkdir(f"{out_file_path}/{system}_poscar")
for idx , dir in enumerate(tqdm(lsqs , desc = "Copying SQS POSCARs")) :
	
	if ".DS" in dir :
		continue
	shutil.copy(
			src = f"{out_file_path}{system}/{dir}/{dir}.vasp" ,
			dst = f"{out_file_path}{system}_poscar/{dir}.vasp"
			)
