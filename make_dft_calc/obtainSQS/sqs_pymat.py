from pymatgen.transformations.advanced_transformations import SQSTransformation
from pymatgen.core.structure import Structure
from pymatgen.core import Element
import warnings
import os, sys


class HiddenPrints:
	def __enter__(self):
		self._original_stdout = sys.stdout
		sys.stdout = open(os.devnull, 'w')
	
	def __exit__(self, exc_type, exc_val, exc_tb):
		sys.stdout.close()
		sys.stdout = self._original_stdout


def obtainSQS(
		inp_file_path: str = '',
		output_path: str = '',
		supercell=None,
		doping_site: str = 'H',
		dopant: str = 'H',
		dopant_percentage: int = 0
):
	"""

	:param struct:
	:param output_path:
	:param supercell:
	:param doping_site:
	:param dopant:
	:param dopant_percentage:
	"""
	if supercell is None:
		supercell = [3, 2, 2]
	struct = Structure.from_file(inp_file_path)
	sqs_transform = SQSTransformation(scaling=supercell, search_time=0.5)
	doping_site = doping_site
	dopant = dopant
	dopant_percentage = dopant_percentage
	for i in struct.sites:
		if i.species.elements[0] == Element(doping_site):
			i.species = {doping_site: 1 - dopant_percentage / 100,
						 dopant: dopant_percentage / 100}
	
	with HiddenPrints():
		warnings.filterwarnings("ignore")
		structure = sqs_transform.apply_transformation(struct)
	
	structure.to(
		f"{output_path}.vasp", fmt="poscar"
	)
