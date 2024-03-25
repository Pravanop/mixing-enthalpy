"""
Creates Phase Diagrams for 325 binaries, 2600 ternaries, 14950 quarternaries, 65780 quinaries
"""
import json
import pickle
from calculateEnthalpy.create_alloy_comp import create_multinary
from pymatgen.core import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram , PDEntry
from tqdm import tqdm

from calculateEnthalpy.thermo_math import gibbs_energy

def make_PD_per_comp(comp , dump_dict , temperature = 1000) :
	"""
	This should be an user-facing function in the future.
	:param comp:
	:param dump_dict:
	:param temperature:
	:return:
	"""
	# comp is "El1-El2-EL3"
	ele_list = comp.split('-')
	n_alloy = len(ele_list)
	
	all_combs = create_multinary(element_list = ele_list , no_comb = list(range(2 , n_alloy + 1)))
	# this is a list of lists
	pd_entry_input = {}
	for idx , n_nary in enumerate(all_combs) :
		idx = idx + 2
		for idx2 , subset_comp in enumerate(n_nary) :
			temp_subset = dump_dict[str(idx)][subset_comp]
			
			# equimolar
			name = Composition(subset_comp.replace('-' , ''))
			pd_entry_input[name] = gibbs_energy(
					temp_subset['mix_enthalpy'] , temp_subset['config_entropy'] ,
					temperature
					) * name.num_atoms
			
			# intermetallics
			if "intermetallic" in temp_subset :
				for idx3 , intermetallic in enumerate(temp_subset['intermetallic']) :
					name = Composition(intermetallic['formula_pretty'])
					pd_entry_input[name] = intermetallic['formation_energy_per_atom'] * name.num_atoms
			
			# elements
			for ele in ele_list :
				name = Composition(ele)
				pd_entry_input[name] = 0 * name.num_atoms
		
		pd_entries_list = []
		for keys , value in pd_entry_input.items() :
			pd_entries_list.append(PDEntry(composition = keys , energy = value))
	
	phase_diagram = PhaseDiagram(pd_entries_list)
	return phase_diagram

def make_phaseDiagram(out_file_name: str = "bokas" , lattice: str = "bcc") :
	"""
	the main function which finishes this step. Creates and stores all the phase diagrams. TODO async multiple
	processes
	:param out_file_name:
	:param lattice:
	"""
	with open(
			"/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/data/output_data/dump_20240217-164122_intermetallic.json" ,
			'r'
			) as f :
		dump_dict = json.load(f)
	
	binary = dump_dict['2']
	ternary = dump_dict['3']
	quartenary = dump_dict['4']
	quinary = dump_dict['5']
	
	phase_diagram_dict = {}
	
	for key , value in tqdm(binary.items() , desc = "Creating Phase Diagrams for all binaries") :
		phase_diagram_dict[key] = make_PD_per_comp(key , dump_dict)
	
	for key , value in tqdm(ternary.items() , desc = "Creating Phase Diagrams for all ternaries") :
		phase_diagram_dict[key] = make_PD_per_comp(key , dump_dict)
	
	for key , value in tqdm(quartenary.items() , desc = "Creating Phase Diagrams for all quartenaries") :
		phase_diagram_dict[key] = make_PD_per_comp(key , dump_dict)
	
	# for key , value in tqdm(quinary.items() , desc = "Creating Phase Diagrams for all quinaries") :
	# 	phase_diagram_dict[key] = make_PD_per_comp(key , dump_dict)
	
	with open(f'data/output_data/PD_{out_file_name}_{lattice}.pickle' , 'w') as f :
		pickle.dump(phase_diagram_dict , f)
	
	print("All phase diagrams have been stored")

# make_phaseDiagram()
