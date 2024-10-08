"""
Creates Phase Diagrams.
Contributor: Joshua Cheng, Pravan Omprakash
"""
import json
import pickle

import numpy as np

from calculateEnthalpy.helper_functions.grid_code import create_multinary
from pymatgen.core import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram , PDEntry
from tqdm import tqdm
import pandas as pd
from calculateEnthalpy.helper_functions.thermo_math import gibbs_energy


def make_PD_per_comp_off_equi(comp, dump_dict, temperature):
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
	# single_energy = load_json_to_dict("./data/input_data/single_energy.json")
	all_combs = create_multinary(element_list=ele_list, no_comb=list(range(2, n_alloy + 1)))
	# this is a list of lists
	pd_entry_input = {}
	for idx, n_nary in enumerate(all_combs):
		idx = idx + 2
		for idx2, subset_comp in enumerate(n_nary):
			temp_subset = dump_dict[str(idx)][subset_comp]
			subset_list = subset_comp.split('-')

			for mol_ratio, value in temp_subset.items():
				if mol_ratio == "intermetallic":
					continue
				mol_ratio = mol_ratio.split('-')
				mol_ratio = np.array(mol_ratio).astype(float)
				mol_ratio = dict(zip(subset_list, mol_ratio))
				mol_ratio = {key: val for key, val in mol_ratio.items() if val != 0.0}
				name = Composition(Composition(mol_ratio).get_integer_formula_and_factor()[0])
				config_entropy = value['config_entropy']
				pd_entry_input[name] = gibbs_energy(
					value['mix_enthalpy'], config_entropy,
					temperature
				) * name.num_atoms

			# intermetallics
			if "intermetallic" in temp_subset:
				for idx3, intermetallic in enumerate(temp_subset['intermetallic']):
					name = Composition(intermetallic['formula_pretty'])
					pd_entry_input[name] = intermetallic['formation_energy_per_atom'] * name.num_atoms

			# elements
			for ele in ele_list:
				name = Composition(ele)
				pd_entry_input[name] = 0 * name.num_atoms

		pd_entries_list = []
		for keys, value in pd_entry_input.items():
			pd_entries_list.append(PDEntry(composition=keys, energy=value))

	phase_diagram = PhaseDiagram(pd_entries_list)
	# decomp, e_hull = phase_diagram.get_decomp_and_e_above_hull(target_pd_entry, allow_negative=True)
	# enthalpy = target_pd_entry.energy_per_atom
	return phase_diagram


def make_phaseDiagram(out_file_name: str = "bokas", lattice: str = "bcc"):
	"""
	the main function which finishes this step. Creates and stores all the phase diagrams. TODO async multiple
	processes
	:param out_file_name:
	:param lattice:
	"""
	with open(
			"/calculateEnthalpy/data/output_data/old_data/dump_20240722-120016_intermetallic.json",
			'r'
	) as f:
		dump_dict = json.load(f)

	binary = dump_dict['2']
	ternary = dump_dict['3']
	quartenary = dump_dict['4']
	quinary = dump_dict['5']

	phase_diagram_dict = {}

	# for key , value in tqdm(binary.items() , desc = "Creating Phase Diagrams for all binaries") :
	# 	phase_diagram_dict[key] = make_PD_per_comp(key , dump_dict)
	for key, value in tqdm(ternary.items(), desc="Creating Phase Diagrams for all ternaries"):
		phase_diagram_dict[key] = make_PD_per_comp(key, dump_dict)

	for key, value in tqdm(quartenary.items(), desc="Creating Phase Diagrams for all quartenaries"):
		phase_diagram_dict[key] = make_PD_per_comp(key, dump_dict)

	for key, value in tqdm(quinary.items(), desc="Creating Phase Diagrams for all quinaries"):
		phase_diagram_dict[key] = make_PD_per_comp(key, dump_dict)
	df = pd.DataFrame({'comp': phase_diagram_dict.keys(), "e_hull": phase_diagram_dict.values()})
	with open(f'data/output_data/PD_{out_file_name}_{lattice}.pickle', 'wb') as f:
		pickle.dump(phase_diagram_dict, f)
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
			"/calculateEnthalpy/data/output_data/old_data/dump_20240806-101354.json",
			'r'
	) as f:
		dump_dict = json.load(f)

	binary = dump_dict['2']
	ternary = dump_dict['3']
	quartenary = dump_dict['4']
	quinary = dump_dict['5']

	phase_diagram_dict = dict()
	temp_dict = dict()
	max_temp = 3200
	# for key , value in tqdm(binary.items() , desc = "Creating Phase Diagrams for all binaries") :
	# 	for temperature in range(1 , max_temp , 200) :
	# 		temp_dict[temperature] = make_PD_per_comp(key , dump_dict , temperature = temperature)
	# 		phase_diagram_dict[key] = temp_dict
	# 	temp_dict = {}
	for key, value in tqdm(ternary.items(), desc="Creating Phase Diagrams for all ternaries"):
		for temperature in range(1, max_temp, 200):
			temp_dict[temperature] = make_PD_per_comp(key, dump_dict, temperature=temperature)
			phase_diagram_dict[key] = temp_dict
		temp_dict = {}
	# for key, value in tqdm(quartenary.items(), desc="Creating Phase Diagrams for all quartenaries"):
	# 	for temperature in range(1, max_temp, 200):
	# 		temp_dict[temperature] = make_PD_per_comp(key, dump_dict, temperature=temperature)
	# 		phase_diagram_dict[key] = temp_dict
	# 	temp_dict = {}
	# for key, value in tqdm(quinary.items(), desc="Creating Phase Diagrams for all quinaries"):
	# 	for temperature in range(1, max_temp, 200):
	# 		temp_dict[temperature] = make_PD_per_comp(key, dump_dict, temperature=temperature)
	# 		phase_diagram_dict[key] = temp_dict
	# 	temp_dict = {}

	with open(f'data/output_data/TD_{out_file_name}_{lattice}.pickle', 'wb') as f:
		pickle.dump(phase_diagram_dict, f)


# make_tempDiagram(out_file_name="john_from_bokas_offequi", lattice="bcc")
