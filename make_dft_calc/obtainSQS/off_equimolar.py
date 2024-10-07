import itertools

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


def obtainSQS_off_equimolar(
		file_path: str,
		system: str,
		elements: list,
		out_folder_path: str
):
	"""

	:param file_path:
	:param system:
	:param elements:
	:param out_folder_path:
	"""
	# structure = Structure.from_file(f"{file_path}")
	
	for doping_percent in [12.5, 25, 50, 75, 87.5]:
		project_name = f"{system}_{doping_percent}"
		
		out_file_path = f"{out_folder_path}{project_name}"
		
		obtainSQS(
			supercell=[3, 2, 1],
			inp_file_path=file_path,
			output_path=out_file_path,
			dopant_percentage=doping_percent,
			dopant=elements[1],
			doping_site=elements[0],
		)


# write_rndstr(
# 		structure = structure ,
# 		output_path = out_file_path ,
# 		doping_percent = doping_percent ,
# 		doping_site = elements[0] ,
# 		doping_ele = elements[1] ,
# 		coord_type = 'frac'
# 		)
#
# write_sqscellout(
# 		supercell = [3 , 2 , 2] ,
# 		output_path = out_file_path
# 		)
#
# sqs_status = runsqs(
# 		corrdump_path = "/Users/pravanomprakash/Documents/Projects/atat/corrdump" ,
# 		mcsqs_path = "/Users/pravanomprakash/Documents/Projects/atat/mcsqs" ,
# 		path_dir = out_file_path ,
# 		pair_interaction = 2.4 ,
# 		triplet_interaction = 3.0
# 		)
#
# sqs2POSCAR(
# 		input_file_path = f"{out_file_path}/bestsqs.out" ,
# 		output_file_path = f"{out_file_path}/{project_name}"
# 		)

element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf']

combs = list(itertools.combinations(element_list, 2))
combs = ['-'.join(sorted(list(i))) for i in combs]
print(combs)
abs_path = "/Users/mcube/Desktop/Projects"
problems = []
for i in combs:
	start, end = i.split('-')
	lattice = "FCC"
	file_path = f"{abs_path}/elements/{lattice}/{start}.vasp"
	out_file_path = f"{abs_path}/make_dft_calc/Outputs_{lattice}/"
	if not os.path.exists(out_file_path):
		os.mkdir(out_file_path)
	if not os.path.exists(f"{out_file_path}{i}"):
		os.mkdir(f"{out_file_path}{i}")
	element_list = [start, end]
	try:
		obtainSQS_off_equimolar(
			file_path=file_path,
			out_folder_path=f"{out_file_path}{i}/",
			system=i,
			elements=element_list,
		)
	except:
		problems.append(i)
		continue

print(problems)
# start = "Ti"
# end = "Hf"
# lattice = "BCC"
# file_path = (f"../Outputs/Outputs_{lattice}/ele_POSCARs/{start}"
# 			 f"_{lattice}.vasp")
# out_file_path = f"../Outputs/Outputs_off_equimolar_{lattice}/"
# system = f"{start}{end}"
# os.mkdir(f"{out_file_path}{system}")
# element_list = [start, end]
# obtainSQS_off_equimolar(
# 	file_path=file_path,
# 	out_folder_path=f"{out_file_path}{system}/",
# 	system=system,
# 	elements=element_list
# )

# lsqs = os.listdir(f"{out_file_path}{system}")
# if not os.path.exists(f"{out_file_path}{system}_poscar/") :
# 	os.mkdir(f"{out_file_path}/{system}_poscar")
# for idx , dir in enumerate(tqdm(lsqs , desc = "Copying SQS POSCARs")) :
#
# 	if ".DS" in dir :
# 		continue
# 	shutil.copy(
# 			src = f"{out_file_path}{system}/{dir}.vasp" ,
# 			dst = f"{out_file_path}{system}_poscar/{dir}.vasp"
# 			)
