from pymatgen.core.structure import Structure
from tqdm import tqdm
from prepareVASPRuns.vasp_run import *

def create_vasprun(
		inp_file_path: str ,
		output_folder_path: str ,
		metadata_path: str = "./metadata/metadata.yaml"
		) :
	"""
	Main Function to create a vasp_run folder. Makes a main folder for ionic steps with name provided by the input
	structure and a "static" subdirectory for the electronic calculations.

	:param inp_file_path: File or folder path. If file path, should be .vasp or .cif file. NO GENERIC NAMES! If
	folder path, then must contain set of .vasp/.cif files only.
	:param output_folder_path: Folder path where output folder is stored
	:param metadata_path: Contains all information needing to create a run. Template is included as default value.
	"""
	metadata = load_yaml_to_dict(metadata_path)
	if identify_file(inp_file_path) :
		loop = [inp_file_path]  # if file
	else :
		loop = os.listdir(inp_file_path)
		loop = [inp_file_path + i for i in loop]  # if folder containing many files
	
	for i in loop:
		
		name = get_filename_without_extension(i)  # name taken from file, hence no generic names.
		named_output_folder_path = f"{output_folder_path}/{name}/"
		create_directory(named_output_folder_path)
		# POSCAR
		species = POSCAR(i , named_output_folder_path)
		
		# KPOINTS
		KPOINTS(metadata , named_output_folder_path , name)
		
		# INCAR
		INCAR(metadata , named_output_folder_path)
		
		# INCAR
		runjob(metadata , named_output_folder_path , name)
		
		# POTCAR
		POTCAR(metadata , species , named_output_folder_path)
		
		# STATIC
		static_folder_path = f"{named_output_folder_path}static/"
		create_directory(static_folder_path)
		
		# Copy POTCAR
		src_file = f"{named_output_folder_path}POTCAR"
		dst_file = f"{static_folder_path}POTCAR"
		copy_file(src_file , dst_file)
		
		# Copy runjob
		src_file = f"{named_output_folder_path}runjob"
		dst_file = f"{static_folder_path}runjob"
		copy_file(src_file , dst_file)
		
		# Increase KPOINT density
		metadata["kpoint"]["grid"] = [element * int(metadata["kpoint"]["static_mult"]) for element in
		                              metadata["kpoint"][
			                              "grid"]]
		KPOINTS(metadata , static_folder_path , name)
		
		# Update INCAR settings. Can add more tags for special calculations.
		metadata["incar"]["update"] = {"ISIF" : "2" , "NSW" : "0" , "EDIFF" : "1E-7"}
		INCAR(metadata , static_folder_path)

