from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram

# Top 5 biggest errors
# Fe-Ta-V, Ta-Ti-W, Ta-Ti-V-W, Cr-Ti-W, Cr-Ta-V

processed_file_path = "/calculateEnthalpy/data/output_data/Pravan_bcc_4/all_lattices_binaries.json"
pD = phaseDiagram(
	processed_file_path=processed_file_path,
grid_size=20)

composition = ["Ti", "Ta", "W"]
equi = True
if equi:
	mol_ratio = [1/len(composition)]*len(composition)


decomp_products, e_above_hull = pD.find_decomp_products(composition=composition,
								  mol_ratio=mol_ratio,
								  temperature=750)
print(decomp_products)