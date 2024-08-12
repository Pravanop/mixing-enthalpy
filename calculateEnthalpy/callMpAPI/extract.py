import os

from pymatgen.core import Element
from tqdm import tqdm

from calculateEnthalpy.callMpAPI.utils import convert_to_conventional


def search_mp_with_EleSG(
		mpr,
		element: str = "C",
		lattice: str = "FCC"
) -> dict:
	"""
    Given the element and space group, it retrieves the structure and the mp_id from MP
    :param element:
    :param lattice:
    :return: dictionary of the form mp_id: conv_s
    """
	spacegroup_lattice_dict = {'BCC': 'Im-3m', 'FCC': 'Fm-3m', 'HCP': 'P6_3/mmc'}
	spacegroup = spacegroup_lattice_dict[lattice]
	docs = mpr.materials.summary.search(
		formula=element,
		spacegroup_symbol=spacegroup,
		fields=["material_id", "structure"]
	)
	try:
		s = docs[0].structure
		mp_id = docs[0].material_id
		conv_s = convert_to_conventional(s)
	except IndexError as e:
		print(e)
		ele = Element(element)
		new_ele = ele.number + 1
		new_ele = Element.from_Z(new_ele)

		docs = mpr.materials.summary.search(
			formula=str(new_ele.name),
			spacegroup_symbol=spacegroup,
			fields=["material_id", "structure"]
		)
		s = docs[0].structure
		mp_id = docs[0].material_id
		conv_s = convert_to_conventional(s)
		for idx, _ in enumerate(conv_s.sites):
			conv_s.replace(species=ele, idx=idx)

	return {element: conv_s}


def search_mp_with_hullE(
		mpr,
		element: str
) -> dict:
	"""
    Given the element , it retrieves the ground state structure and the mp_id from MP
    :param element:
    :param lattice:
    :return: dictionary of the form mp_id: conv_s
    """
	docs = mpr.materials.summary.search(
		formula=element,
		energy_above_hull=0.0,
		fields=["material_id", "structure", "spacegroup"]
	)
	s = docs[0].structure
	mp_id = docs[0].material_id
	conv_s = convert_to_conventional(s)
	return {element: conv_s}


"""
Step 3 is to create a for loop to run across the entire list.
"""


def batch_extract_EleSG(
		mpr,
		element_list: list[str],
		lattice: str = "FCC"
) -> dict:
	"""
    
    :param element_list: list of element symbols
    :param lattice: the lattice system that is being extracted
    :return: dictionary of (mp_id: structure) pairs
    """
	results_dict = {}
	for idx, ele in tqdm(enumerate(element_list), desc="Extracting structures for element list"):
		ele_dict = search_mp_with_EleSG(
			mpr=mpr,
			element=ele,
			lattice=lattice
		)
		results_dict.update(ele_dict)

	return results_dict


def batch_extract_hullE(
		mpr,
		element_list: list[str],
) -> dict:
	"""

    :param element_list: list of element symbols
    :param lattice: the lattice system that is being extracted
    :return: dictionary of (mp_id: structure) pairs
    """
	results_dict = {}
	for idx, ele in tqdm(enumerate(element_list), desc="Extracting structures for element list"):
		ele_dict = search_mp_with_hullE(
			mpr=mpr,
			element=ele,
		)
		results_dict.update(ele_dict)

	return results_dict


"""
Step 4 is to take this dict and store the POSCAR files in a specific location
"""


def eleList_to_POSCAR(
		mpr,
		element_list: list[str],
		lattice: str,
		extract_style: str,
		out_file_path: str
) -> dict:
	"""
    
    :param element_list: list of element symbols
    :param lattice: the lattice system that is being extracted
    :param out_file_path: folder path to store POSCAR files
    :return: element_dict if needed for further processing or checking
    """
	ele_dict = {}
	if extract_style == "sym":
		ele_dict = batch_extract_EleSG(mpr=mpr, element_list=element_list, lattice=lattice)
	if extract_style == "ground":
		ele_dict = batch_extract_hullE(mpr=mpr, element_list=element_list)

	if not os.path.exists(out_file_path):
		os.mkdir(out_file_path)

	for key, value in tqdm(ele_dict.items(), desc="Storing POSCAR files of extracted elements"):
		file_path = f"{out_file_path}{str(key)}_{lattice}.vasp"
		value.to(fmt='poscar', filename=file_path)

	return ele_dict
