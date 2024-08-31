from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram

processedbinary__file_path = "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/new_phase_diagram/bokas_omegas_processed.json"
end_member_path = "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/new_phase_diagram/bokas_end_members_dict.json"
pD = phaseDiagram(
	processed_binary_file_path=processedbinary__file_path,
	end_member_file_path=end_member_path,
grid_size=30)

composition = ['W', 'Ti']
equi = True
if equi:
	mol_ratio = [1/len(composition)]*len(composition)

for lattice in ['BCC', 'FCC', 'HCP']:
	misc_T = pD.find_misc_temperature(composition=composition, lattice=lattice, mol_ratio=mol_ratio, flag='off_equi', correction=False)
	print(lattice, misc_T)
