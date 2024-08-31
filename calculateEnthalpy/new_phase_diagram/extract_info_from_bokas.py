import json

import numpy as np

with open("./bokas_omegas.json", 'r') as f:
	data = json.load(f)

binary_omegas = data['omegas']
end_members = data['elements']

# print(binary_omegas.keys(), end_members.keys())

def calc_delta_H(end_members_dict):
	lattices = list(end_members_dict.keys())
	members_dict = {}
	for lattice in lattices:
		for member in end_members_dict[lattice]:
			if member not in members_dict:
				members_dict[member] = {lattice: end_members_dict[lattice][member]}
			else:
				members_dict[member].update({lattice: end_members_dict[lattice][member]})

	for key, value in members_dict.items():
		lattice_energies = np.array(list(value.values()))
		min_lattice_energy = min(lattice_energies)
		lattice_energies = np.round(lattice_energies - min_lattice_energy,4) # in eV/atom
		members_dict[key] = dict(zip(value.keys(), lattice_energies))
	return members_dict

members_dict = calc_delta_H(end_members)

#storing the processed end_member_file
with open('./bokas_end_members_dict.json', 'w') as f:
	json.dump(members_dict, f, indent=4, ensure_ascii=True)

#processing the omegas too.

def process_omegas(omegas_dict):
	lattices = list(omegas_dict.keys())
	omegas_dict_processed = {}
	for lattice in lattices:
		for member in omegas_dict[lattice]:
			if member not in omegas_dict_processed:
				omegas_dict_processed[member] = {lattice: np.round(omegas_dict[lattice][member],4)}
			else:
				omegas_dict_processed[member].update({lattice: np.round(omegas_dict[lattice][member],4)})

	return omegas_dict_processed

omegas = process_omegas(binary_omegas)
with open('./bokas_omegas_processed.json', 'w') as f:
	json.dump(omegas, f, indent=4, ensure_ascii=True)