import json
import pickle
from create_alloy_comp import create_multinary
from pymatgen.core import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram , PDEntry
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt

from thermo_math import gibbs_energy

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
			if subset_comp == comp:
				target_pd_entry = PDEntry(name, gibbs_energy(
					temp_subset['mix_enthalpy'], temp_subset['config_entropy'],
					temperature
				) * name.num_atoms)
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
	e_hull = phase_diagram.get_e_above_hull(target_pd_entry, allow_negative=True)
	return e_hull

def make_phaseDiagram(out_file_name: str = "bokas" , lattice: str = "bcc") :
	"""
	the main function which finishes this step. Creates and stores all the phase diagrams. TODO async multiple
	processes
	:param out_file_name:
	:param lattice:
	"""
	with open(
			"data/output_data/dump_20240322-125557_intermetallic.json" ,
			'r'
			) as f :
		dump_dict = json.load(f)
	
	binary = dump_dict['2']
	ternary = dump_dict['3']
	quartenary = dump_dict['4']
	quinary = dump_dict['5']
	
	phase_diagram_dict = {}
	
	# for key , value in tqdm(binary.items() , desc = "Creating Phase Diagrams for all binaries") :
	# 	phase_diagram_dict[key] = make_PD_per_comp(key , dump_dict)
	for key , value in tqdm(ternary.items() , desc = "Creating Phase Diagrams for all ternaries") :
		phase_diagram_dict[key] = make_PD_per_comp(key , dump_dict)

	for key , value in tqdm(quartenary.items() , desc = "Creating Phase Diagrams for all quartenaries") :
		phase_diagram_dict[key] = make_PD_per_comp(key , dump_dict)

	for key , value in tqdm(quinary.items() , desc = "Creating Phase Diagrams for all quinaries") :
		phase_diagram_dict[key] = make_PD_per_comp(key , dump_dict)
	df = pd.DataFrame({'comp': phase_diagram_dict.keys(), "e_hull": phase_diagram_dict.values()})
	with open(f'data/output_data/PD_{out_file_name}_{lattice}.pickle' , 'wb') as f :
		pickle.dump(phase_diagram_dict , f)
	df.to_excel('test.xlsx', index=False)
	print("All phase diagrams have been stored")


def make_tempDiagram(out_file_name: str = "bokas", lattice: str = "bcc"):
	"""
	the main function which finishes this step. Creates and stores all the phase diagrams. TODO async multiple
	processes
	:param out_file_name:
	:param lattice:
	"""
	with open(
			"data/output_data/dump_20240322-125557_intermetallic.json",
			'r'
	) as f:
		dump_dict = json.load(f)

	binary = dump_dict['2']
	ternary = dump_dict['3']
	quartenary = dump_dict['4']
	quinary = dump_dict['5']

	phase_diagram_dict = dict()
	temp_dict = dict()
	for key, value in tqdm(ternary.items(), desc="Creating Phase Diagrams for all ternaries"):
		for temperature in range(0, 2200, 200):
			temp_dict[temperature] = make_PD_per_comp(key, dump_dict, temperature=temperature)
			phase_diagram_dict[key] = temp_dict
		temp_dict = {}
	for key, value in tqdm(quartenary.items(), desc="Creating Phase Diagrams for all quartenaries"):
		for temperature in range(0, 2200, 200):
			temp_dict[temperature] = make_PD_per_comp(key, dump_dict, temperature=temperature)
			phase_diagram_dict[key] = temp_dict
		temp_dict = {}
	for key, value in tqdm(quinary.items(), desc="Creating Phase Diagrams for all quinaries"):
		for temperature in range(0, 2200, 200):
			temp_dict[temperature] = make_PD_per_comp(key, dump_dict, temperature=temperature)
			phase_diagram_dict[key] = temp_dict
		temp_dict = {}
	with open(f'data/output_data/TD_{out_file_name}_{lattice}.json' , 'w') as f :
		json.dump(phase_diagram_dict , f)

# make_tempDiagram()

def plot_tempDiagram(comp):
	with open(
			"data/output_data/TD_bokas_bcc.json",
			'r'
	) as f:
		dump_dict = json.load(f)
	data = dump_dict[comp]
	for key, value in data.items():
		data[key] = value * 1000
	temperatures = data.keys()
	energy_above_hull = data.values()
	print(temperatures)
	print(energy_above_hull)
	# Plotting
	plt.figure(figsize=(8, 6))
	plt.plot(temperatures, energy_above_hull, marker='o', linestyle='-')

	# Customizing plot appearance
	plt.xlabel('Temperature (K)')
	plt.ylabel('Energy above Hull (meV)')
	plt.title(comp + " Stability Diagram")
	plt.grid(True)

	# Save or show the plot
	plt.tight_layout()  # Adjust layout to prevent clipping of labels
	plt.savefig('phase_diagram.png', dpi=300)  # Save plot as image
	plt.show()

plot_tempDiagram("Cr-Fe-V")