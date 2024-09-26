import os
from mp_api.client import MPRester
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

lattice = 'HCP'
os.mkdir(f"../elements/{lattice}/")
element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf', 'Fe']
for i in element_list:
	mpr = MPRester(api_key = "u1TjwfwfTnpF8IolXF9PBY9RT9YauL84" , mute_progress_bars = True)
	ans = search_mp_with_EleSG(mpr , element = i , lattice = lattice)
	ans[i].to(f"../elements/{lattice}/{i}.vasp" , fmt = "poscar")
