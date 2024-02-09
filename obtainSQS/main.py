import numpy as np
from pymatgen.core.structure import Structure

def load_structure(path: str = ""):
	return Structure.from_file(path)

path = "/Users/pravanomprakash/Library/CloudStorage/Box-Box/LaHfO2/FE_02_111_rel.vasp"

def get_parameters(path):
	s = load_structure(path = path)
	positions = s.frac_coords
	lattice = s.lattice.matrix
	elements = s.sites
	return lattice, positions, elements

doping_percent = 12.5
doping_ele = "La"
doping_site = "Hf"

def write_rndstr(path, supercell, output_path):
	params = get_parameters(path)
	#lattice lines
	lines = ""
	lines += "1 1 1 90 90 90\n"
	
	#supercell lines
	lines += f"{supercell} 0.0 0.0\n"
	lines += f"0.0 {supercell} 0.0\n"
	lines += f"0.0 0.0 {supercell}\n"
	
	#coordinates lines
	for i in params[2]:
		if str(i.specie) == doping_site:
			lines += (f"{np.round(i.x,3)} {np.round(i.y,3)} {np.round(i.z,3)} {i.specie}={1 - doping_percent/100}"
			          f",{doping_ele}"
			          f"={doping_percent/100}\n")
		else:
			lines += f"{np.round(i.x,3)} {np.round(i.y,3)} {np.round(i.z,3)} {i.specie}\n"
	
	with open(f'{output_path}/rndstr.in' , 'w') as file :
		file.write(lines)

def write_conc(path, supercell, output_path):
	params = get_parameters(path)
	# lattice lines
	lines = ""
	for i in params[0] :
		lines += f"{np.round(i[0] , 2)} {np.round(i[1] , 2)} {np.round(i[2] , 2)}\n"
	
	# supercell lines
	lines += f"{supercell} 0.0 0.0\n"
	lines += f"0.0 {supercell} 0.0\n"
	lines += f"0.0 0.0 {supercell}\n"
	
	# coordinates lines
	for i in params[2] :
		lines += f"{np.round(i.x , 3)} {np.round(i.y , 3)} {np.round(i.z , 3)} {i.specie}\n"
	
	with open(f'{output_path}/conc.in' , 'w') as file :
		file.write(lines)
		
write_rndstr(path, 2, output_path = "/Users/pravanomprakash/Downloads/lahfo2_12")
write_conc(path, 2, output_path = "/Users/pravanomprakash/Downloads/lahfo2_12")


