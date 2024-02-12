from pymatgen.core.structure import Structure
from prepareVASPRuns.file_utils import *
def POSCAR(inp_file_path, named_output_folder_path):
	s = Structure.from_file(inp_file_path)
	s.to(f"{named_output_folder_path}POSCAR" , fmt = "poscar")
	species = get_element_symbols(s)
	return species

def KPOINTS(metadata, named_output_folder_path, name):
	grid = metadata['kpoint']['grid']
	# print(type(grid))
	lines = ""
	lines += f"Comment: {metadata['kpoint']['type']} grid for {name}\n"
	lines += "0           Automatic generation of k-points\n"
	lines += f"{metadata['kpoint']['type']}\n"
	lines += f"{grid[0]} {grid[1]} {grid[0]}\n"
	kpoints_file_path = f"{named_output_folder_path}KPOINTS"
	string_to_file(lines , kpoints_file_path)
	
	
def INCAR(metadata, named_output_folder_path):
	incar_dict = load_json_to_dict(metadata["incar"]["path"])
	if metadata["incar"]["update"] != 'None' :
		incar_dict.update(metadata["incar"]["update"])
	incar_file = ""
	for key , value in incar_dict.items() :
		incar_file += f"{key} = {value}\n"
	incar_file_path = f"{named_output_folder_path}INCAR"
	string_to_file(incar_file , incar_file_path)
	

def runjob(metadata, named_output_folder_path, name):
	runjob_meta = metadata["runjob"]
	runjob_file_path = f"{named_output_folder_path}runjob"
	if runjob_meta["supercomputer"] == "expanse" :
		runjob_lines = ""
		runjob_lines += f"#!/bin/bash\n"
		runjob_lines += f"#SBATCH --job-name={name}\n"
		runjob_lines += f"#SBATCH --partition=\{runjob_meta['comp_name']}\n"
		runjob_lines += f"#SBATCH --time 0{runjob_meta['hours']}:00:00\n"
		runjob_lines += f"#SBATCH --output=vasp.out\n"
		runjob_lines += f"#SBATCH --nodes={runjob_meta['nodes']}\n"
		runjob_lines += f"#SBATCH --ntasks-per-node={runjob_meta['cores']}\n"
		runjob_lines += f"#SBATCH --account=tg-dmr160007\n"
		runjob_lines += f"#SBATCH --mem=0\n"
		runjob_lines += f"module purge\nmodule load cpu/0.15.4\n"
		runjob_lines += f"module load  intel/19.1.1.217\nmodule load slurm\n"
		runjob_lines += f"module load intel-mkl/2019.1.144\n"
		runjob_lines += f"module load intel-mpi/2019.8.254\n"
		runjob_lines += f"mpirun $VASP >> stdout\n"
		string_to_file(runjob_lines , runjob_file_path)
	
	elif runjob_meta["supercomputer"] == "stampede" :
		runjob_lines = ""
		runjob_lines += f"#!/bin/bash\n"
		runjob_lines += f"#SBATCH -J {name}\n"
		runjob_lines += f"#SBATCH -o status.o%j \n"
		runjob_lines += f"#SBATCH -e status.e%j \n"
		runjob_lines += f"#SBATCH -p {runjob_meta['comp_name']}\n"
		runjob_lines += f"#SBATCH -N {runjob_meta['nodes']}\n"
		runjob_lines += f"#SBATCH -n {runjob_meta['cores']}\n"
		runjob_lines += f"#SBATCH -t 0{runjob_meta['hours']}:00:00\n"
		runjob_lines += f"#SBATCH -A TG-DMR160007\n"
		runjob_lines += f"export FORT_BUFFERED=TRUE\n"
		runjob_lines += f"module load intel/18.0.2\n"
		runjob_lines += f"ibrun $VASP >> stdout\n"
		string_to_file(runjob_lines , runjob_file_path)

def POTCAR(metadata, species, named_output_folder_path):
	if metadata['potcar']['functional'] == 'LDA' :
		potcar_folder_path = f"{metadata['potcar']['path']}/potpaw"
	
	else :
		potcar_folder_path = f"{metadata['potcar']['path']}/potpaw_{metadata['potcar']['functional']}"
	
	potcar_dict = load_yaml_to_dict(metadata["potcar"]["base_path"])['POTCAR']

	elem_paths = []
	for elem in species :
		elem_paths.append(f"{potcar_folder_path}/{potcar_dict[elem]}/POTCAR")
	
	potcar_file_path = f"{named_output_folder_path}POTCAR"
	concatenate_files(elem_paths , potcar_file_path)
	