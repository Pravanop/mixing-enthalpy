from pymatgen.transformations.advanced_transformations import SQSTransformation
from pymatgen.core.structure import Structure
from pymatgen.core import Element

def obtainSQS(
		input_path: str = '' ,
		output_path: str = '' ,
		supercell = None ,
		doping_site: str = 'H' ,
		dopant: str = 'H' ,
		dopant_percentage: int = 0
		) :
	if supercell is None :
		supercell = [3 , 2 , 2]
	sqs_transform = SQSTransformation(scaling = supercell , search_time = 1)
	s = Structure.from_file(
		filename = input_path
		)
	doping_site = doping_site
	dopant = dopant
	dopant_percentage = dopant_percentage
	for i in s.sites :
		if i.species.elements[0] == Element(doping_site) :
			i.species = {doping_site : 1 - dopant_percentage / 100 , dopant : dopant_percentage / 100}
	
	structure = sqs_transform.apply_transformation(s)
	structure.to(
		output_path , fmt = "poscar"
		)
