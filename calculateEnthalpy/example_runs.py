import json

from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram
from calculateEnthalpy.helper_functions.thermo_math import thermoMaths

processed_file_path = "/calculateEnthalpy/data/output_data/Pravan_bcc_4/all_lattices_binaries.json"
pD = phaseDiagram(
	processed_file_path=processed_file_path,
grid_size=30)
composition = ['Cr', 'W']
equi = True
if equi:
	mol_ratio = [1/len(composition)]*len(composition)

# #Find Miscibility temperature
misc_T = pD.find_misc_temperature(composition=composition, mol_ratio=mol_ratio, flag='equi')
print(misc_T)
if misc_T is None:
	misc_T = 3200

if misc_T != 0:
	decomp_products, e_above_hull = pD.find_decomp_products(composition=composition,
								  mol_ratio=mol_ratio,
								  temperature=2200)

	print(f'E_above_hull: at {2200}K ', e_above_hull, " eV/atom")
	print(decomp_products)
#
#Find mixing enthalpy of a multinary
tm = thermoMaths()
composition_dict =  dict(zip(composition, mol_ratio))
composition_dict = {key: val for key, val in composition_dict.items() if val != 0.0}
with open(processed_file_path, "r") as f:
	data = json.load(f)
print(tm.calc_mutinary_multilattice_mix_Enthalpy(composition_dict, data))

#Find mixing enthalpies of all subsets in a composition
print(tm.find_subset_enthalpies(composition, data))



