import json
import pickle

import numpy as np

from create_alloy_comp import create_multinary
from pymatgen.core import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram , PDEntry
from tqdm import tqdm
import pandas as pd
from prepareVASPRuns.file_utils import load_json_to_dict
from thermo_math import gibbs_energy , vibrational_energy

def make_PD_per_comp(comp , dump_dict , temperature = 1000) :
	"""
	This should be a user-facing function in the future.
	:param comp:
	:param dump_dict:
	:param temperature:
	:return:
	"""
	# comp is "El1-El2-EL3"
	ele_list = comp.split('-')
	n_alloy = len(ele_list)
	single_energy = load_json_to_dict("./data/input_data/single_energy.json")
	all_combs = create_multinary(element_list = ele_list , no_comb = list(range(2 , n_alloy + 1)))
	# this is a list of lists
	pd_entry_input = {}
	for idx , n_nary in enumerate(all_combs) :
		idx = idx + 2
		for idx2 , subset_comp in enumerate(n_nary) :
			temp_subset = dump_dict[str(idx)][subset_comp]
			# print(subset_comp)
			# equimolar
			name = Composition(subset_comp.replace('-' , ''))
			subset_comp_list = subset_comp.split('-')
			vib_energy = []
			for e in subset_comp_list :
				e_debye = single_energy[e]["debye"]
				temp_vib = vibrational_energy(
					debye_temp = e_debye ,
					temp = temperature
					)
				vib_energy.append(temp_vib)
			avg_vib_energy = np.mean(np.array(vib_energy))
			pd_entry_input[name] = gibbs_energy(
					temp_subset['mix_enthalpy'] , temp_subset['config_entropy'] ,
					temperature
					) * name.num_atoms
			if subset_comp == comp :
				target_pd_entry = PDEntry(
					name , gibbs_energy(
							temp_subset['mix_enthalpy'] , temp_subset['config_entropy'] ,
							temperature
							) * name.num_atoms
					)
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
	decomp , e_hull = phase_diagram.get_decomp_and_e_above_hull(target_pd_entry , allow_negative = True)
	enthalpy = target_pd_entry.energy_per_atom
	return e_hull , decomp , enthalpy

def make_phaseDiagram(out_file_name: str = "bokas" , lattice: str = "bcc") :
	"""
	the main function which finishes this step. Creates and stores all the phase diagrams. TODO async multiple
	processes
	:param out_file_name:
	:param lattice:
	"""
	with open(
			"data/output_data/dump_20240401-112223_intermetallic.json" ,
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
	df = pd.DataFrame({'comp' : phase_diagram_dict.keys() , "e_hull" : phase_diagram_dict.values()})
	with open(f'data/output_data/PD_{out_file_name}_{lattice}.pickle' , 'wb') as f :
		pickle.dump(phase_diagram_dict , f)
	df.to_excel('test.xlsx' , index = False)
	print("All phase diagrams have been stored")

def make_tempDiagram(out_file_name: str = "bokas" , lattice: str = "bcc") :
	"""
	the main function which finishes this step. Creates and stores all the phase diagrams. TODO async multiple
	processes
	:param out_file_name:
	:param lattice:
	"""
	with open(
			"data/output_data/dump_20240401-112223_intermetallic.json" ,
			'r'
			) as f :
		dump_dict = json.load(f)
	
	binary = dump_dict['2']
	ternary = dump_dict['3']
	quartenary = dump_dict['4']
	quinary = dump_dict['5']
	
	phase_diagram_dict = dict()
	temp_dict = dict()
	max_temp = 3200
	for key , value in tqdm(binary.items() , desc = "Creating Phase Diagrams for all binaries") :
		for temperature in range(1 , max_temp , 200) :
			temp_dict[temperature] = make_PD_per_comp(key , dump_dict , temperature = temperature)
			phase_diagram_dict[key] = temp_dict
		temp_dict = {}
	for key , value in tqdm(ternary.items() , desc = "Creating Phase Diagrams for all ternaries") :
		for temperature in range(1 , max_temp , 200) :
			temp_dict[temperature] = make_PD_per_comp(key , dump_dict , temperature = temperature)
			phase_diagram_dict[key] = temp_dict
		temp_dict = {}
	for key , value in tqdm(quartenary.items() , desc = "Creating Phase Diagrams for all quartenaries") :
		for temperature in range(1 , max_temp , 200) :
			temp_dict[temperature] = make_PD_per_comp(key , dump_dict , temperature = temperature)
			phase_diagram_dict[key] = temp_dict
		temp_dict = {}
	for key , value in tqdm(quinary.items() , desc = "Creating Phase Diagrams for all quinaries") :
		for temperature in range(1 , max_temp , 200) :
			temp_dict[temperature] = make_PD_per_comp(key , dump_dict , temperature = temperature)
			phase_diagram_dict[key] = temp_dict
		temp_dict = {}

	with open(f'data/output_data/TD_{out_file_name}_{lattice}.pickle' , 'wb') as f :
		pickle.dump(phase_diagram_dict , f)

make_tempDiagram()
