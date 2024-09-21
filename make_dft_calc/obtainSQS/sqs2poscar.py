"""
Parser to convert bestsqs.out from mcsqs command in ATAT to POSCAR.
Contributer: Pravan Omprakash
"""

import numpy as np
from pymatgen.core.structure import Structure

def sqs2POSCAR(
		input_file_path: str ,
		output_file_path: str
		) :
	with open(input_file_path , 'r') as f :
		sqs = f.read()
	
	sqs_as_list = sqs.splitlines()
	
	lattice = sqs_as_list[:3]
	lattice = np.array([x.split(' ') for x in lattice]).astype(float)
	supercell = sqs_as_list[3 :6]
	supercell = np.array([x.split(' ') for x in supercell]).astype(float)
	
	lattice_final = supercell * lattice
	lattice_final = lattice_final
	
	positions = sqs_as_list[6 :]
	positions = np.array([x.split(' ') for x in positions])
	elements = positions[: , -1]
	cart_coords = np.array(positions[: , :-1]).astype(float)
	cart_coords_final = np.matmul(cart_coords , lattice)
	
	S = Structure(
		lattice = lattice_final ,
		coords_are_cartesian = True ,
		species = elements ,
		coords = cart_coords_final
		)
	S.to_file(f'{output_file_path}.vasp', fmt = 'poscar')
	
	return True
