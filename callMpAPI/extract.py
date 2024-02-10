from mp_api.client import MPRester
import os
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure
from tqdm import tqdm

from utils import *

"""Goal of this script is to take in a set of elements, a specific spacegroup and return the structure,
mpid and other useful information"""

"""First let us define a few helper functions. We will move these functions to an utils file when finished with
project."""




"""
Step 1 is to define the inputs that are needed from user
"""

element_list = []
spacegroup = ""
xc_functional = "GGA_GGA+U"
spacegroup_lattice_dict = {'BCC': 'Im-3m', 'FCC': 'Fm-3m', 'HCP': 'P6_3/mmc'}
mpr = MPRester(api_key=getAPIKey())

"""
Step 2 is to make a simple function for extracting just structure and material id for a single element.
"""
def search_mp_with_EleSG(element: str = "C",
                         lattice: str = "FCC") -> dict:
    """
    Given the element and space group, it retrieves the structure and the mp_id from MP
    :param element:
    :param lattice:
    :return: dictionary of the form mp_id: conv_s
    """
    spacegroup = spacegroup_lattice_dict[lattice]
    docs = mpr.materials.summary.search(
        formula = element ,
        spacegroup_symbol = spacegroup ,
        fields = ["material_id", "structure"]
        )
    s = docs[0].structure
    mp_id = docs[0].material_id
    conv_s = convert_to_conventional(s)
    return {mp_id : conv_s}

"""
Step 3 is to create a for loop to run across the entire list.
"""

def batch_extract_EleSG(element_list : list[str] = ["C"],
                        lattice: str = "FCC") -> dict:
    """
    
    :param element_list: list of element symbols
    :param lattice: the lattice system that is being extracted
    :return: dictionary of (mp_id: structure) pairs
    """
    results_dict = {}
    for idx, ele in tqdm(enumerate(element_list), desc = "Extracting structures for element list"):
        
        ele_dict = search_mp_with_EleSG(element = ele,
                                   lattice = lattice)
        results_dict.update(ele_dict)
    
    return results_dict

"""
Step 4 is to take this dict and store the POSCAR files in a specific location
"""

def eleList_to_POSCAR(element_list : list[str] = ["C"],
                      lattice: str = "FCC",
                      out_file_path:str = "./poscar_files_step1/") -> dict:
    """
    
    :param element_list: list of element symbols
    :param lattice: the lattice system that is being extracted
    :param out_file_path: folder path to store POSCAR files
    :return: element_dict if needed for further processing or checking
    """
    ele_dict = batch_extract_EleSG(element_list = element_list, lattice = "FCC")
    if not os.path.exists(out_file_path) :
        os.mkdir(out_file_path)
    
    for key, value in tqdm(ele_dict.items(), desc = "Storing POSCAR files of extracted elements"):
        file_path =  f"{out_file_path}{str(key)}_{lattice}.vasp"
        value.to(fmt='poscar', filename=file_path)
        
    
 ## just for testing
if __name__ == '__main__':
    els = ['Al' , 'Co' , 'Cr' , 'Cu' , 'Fe' , 'Hf' , 'Mn' , 'Mo' ,
               'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr', 'V', 'Mg', 'Re',
               'Os', 'Rh', 'Ir', 'Pd', 'Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Ru']
    eleList_to_POSCAR(element_list = els,
                      lattice = "FCC")

##Check list:
#TODO Error handling
#TODO Get ground state phases instead of specific spacegroup
#TODO suppress mp output for cleaner print statements


##old code
# # temporary
# a = getAPIKey()
# mpr = MPRester(api_key=a)
# with MPRester(api_key=a) as mpr:
#     thermo_docs = mpr.materials.search(
#         material_ids=["mp-8635"], fields=["material_id", "initial_structures", "symmetry"]
#     )
#
# # thermo_types=["GGA_GGA+U"]
#
# # takes list of elements and runs structure_search to retrieve necessary info into dict
# def list_to_dict(element_list, spacegroup:str, thermo='GGA_GGA+U'):
#     out = {}
#     for element in element_list:
#         mp_id = structure_search(element, spacegroup)
#         if mp_id == -1:
#             print("No Structure Found")
#         else:
#             thermo_docs = mpr.thermo.search(
#                 material_ids=[mp_id], thermo_types=[thermo], fields=["energy_above_hull"]
#             )
#             out[element] = thermo_docs[0].energy_above_hull
#     return out
#
# # searches for summary based on element
# def structure_search(element: str, spacegroup: str = ""):
#     # if spacegroup isn't given (should be)
#     if not spacegroup:
#         docs = mpr.summary.search(formula=element, fields=["material_id", "symmetry"])
#         return docs
#     if isinstance(spacegroup, str):
#         struct_alias = {'bcc': 'Im-3m', 'fcc': 'Fm-3m', 'hcp': 'P6_3/mmc'}
#         if spacegroup in struct_alias.keys():
#             spacegroup = struct_alias[spacegroup]
#         docs = mpr.summary.search(formula=element, spacegroup_symbol=spacegroup, fields=["structure", "symmetry", "material_id"])
#         print(docs)
#         if not docs:
#             return -1
#         structure = docs[0].structure
#         directory = "POSCAR_FILES/"
#         if not os.path.exists(directory):
#             os.mkdir(directory)
#         structure.to(fmt='poscar', filename=directory+element+"_"+spacegroup)
#         return docs[0].material_id
#     return -1
#
# # if __name__ == '__main__':
# #     els = ['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn', 'Mo',
# #            'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr', 'V', 'Mg', 'Re',
# #            'Os', 'Rh', 'Ir', 'Pd', 'Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Ru']
# #     print(list_to_dict(els, 'fcc'))
# print(thermo_docs)