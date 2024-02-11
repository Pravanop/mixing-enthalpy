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