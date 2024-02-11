import numpy as np
import pymatgen.core.structure
from pymatgen.core.structure import Structure

"""
Get SQS input files from pymatgen structure as required for commands "corrdump" and "mcsqs" in ATAT.
Contributor : Pravan Omprakash
"""

def load_structure(path: str = "") :
	"""
	
	:param path:
	:return:
	"""
	return Structure.from_file(path)

def get_parameters(s : pymatgen.core.structure.Structure) :
	"""
	
	:param s:
	:return:
	"""
	positions = s.frac_coords
	lattice = s.lattice.matrix
	elements = s.sites
	return lattice , positions , elements

def write_rndstr(
		structure: pymatgen.core.structure.Structure ,
		output_path: str ,
		doping_percent: float ,
		doping_ele: str ,
		doping_site: str
		) :
	"""
	
	:param structure:
	:param supercell:
	:param output_path:
	:param doping_percent:
	:param doping_ele:
	:param doping_site:
	:return:
	"""
	
	params = get_parameters(structure)
	
	# lattice lines
	lines = ""
	for i in params[0] :
		lines += f"{np.round(i[0] , 2)} {np.round(i[1] , 2)} {np.round(i[2] , 2)}\n"
	
	# supercell lines
	lines += f"1 0.0 0.0\n"
	lines += f"0.0 1 0.0\n"
	lines += f"0.0 0.0 1\n"
	
	# coordinates lines
	for i in params[2] :
		if str(i.specie) == doping_site :
			lines += (
					f"{np.round(i.x , 3)} {np.round(i.y , 3)} {np.round(i.z , 3)} {i.specie}="
					f"{1 - doping_percent / 100}"
					f",{doping_ele}"
					f"={doping_percent / 100}\n")
		else :
			lines += f"{np.round(i.x , 3)} {np.round(i.y , 3)} {np.round(i.z , 3)} {i.specie}\n"
	
	with open(f'{output_path}/rndstr.in' , 'w') as file :
		file.write(lines)


def write_sqscellout(
		supercell: list[int , int , int] ,
		output_path: str ,
		) :
	lines = ""
	lines += "1\n"
	lines += f"{supercell[0]} 0 0\n"
	lines += f"0 {supercell[1]} 0\n"
	lines += f"0 0 {supercell[2]}\n"
	
	with open(f'{output_path}/sqscell.out' , 'w') as file :
		file.write(lines)
