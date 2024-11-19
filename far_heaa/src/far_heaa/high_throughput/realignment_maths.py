import numpy as np
import itertools
from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations

def compute_plane_4d(*points):
	n = len(points[0])  # Dimensionality of the space
	if len(points) != n - 1:
		raise ValueError("Number of points must be n-1 for an n-dimensional hyperplane.")
	
	# Compute vectors spanning the hyperplane
	base_point = points[0]
	vectors = [point - base_point for point in points[1:]]
	
	# Stack the vectors into a matrix
	matrix = np.stack(vectors)
	
	# Find the null space of the matrix (normal vector to the hyperplane)
	u, s, vh = np.linalg.svd(matrix, full_matrices=True)
	normal_vector = vh[-1]  # The last row of V^H is the null space
	
	# Compute the constant term d for the hyperplane equation
	d = -np.dot(normal_vector, base_point)
	
	return normal_vector, d


def compute_energy_from_plane(normal_vector, d, mol_fractions):
	"""
	Compute the energy (E) from a plane equation given mol fractions.

	Parameters:
		normal_vector (numpy array): The normal vector of the plane [n_x, n_y, ..., n_E].
		d (float): The constant term in the hyperplane equation.
		mol_fractions (numpy array): The mol fractions [x, y, z, ...].

	Returns:
		float: The energy (E) calculated from the plane.
	"""
	# Extract the coefficients for the mol fractions (all except the last)
	mol_fraction_coeffs = normal_vector[:-1]
	
	# Coefficient for energy
	energy_coeff = normal_vector[-1]
	
	# Calculate energy using the hyperplane equation
	numerator = - (np.dot(mol_fraction_coeffs, mol_fractions) + d)
	energy = numerator / energy_coeff
	
	return energy


def get_realignment_energy(element_list, lattice, data):
	n = len(element_list) - 1
	mol_n_list = []
	for i in range(n):
		temp = np.ones(n + 1)
		temp[i] = 0
		mol_n_list.append(temp / n)
	
	mol_n_list = np.array(mol_n_list)
	
	binary_pairs = list(MultinaryCombinations.create_multinary(element_list=element_list, no_comb=[2]).values())[0]
	H_full = np.array([data[i][lattice] for i in binary_pairs])
	realign_vectors = []
	for i in mol_n_list:
		mol_comb = np.array(list(itertools.combinations(i, r=2)))
		comp_product = mol_comb[:, 0] * mol_comb[:, 1]
		H_c_all = np.sum(comp_product * H_full)
		vector = np.append(i[:-1], H_c_all)
		realign_vectors.append(vector)
	
	total_vector = np.array([1 / (n + 1)] * (n + 1))
	mol_comb = np.array(list(itertools.combinations(total_vector, r=2)))
	comp_product = mol_comb[:, 0] * mol_comb[:, 1]
	H_c_all = np.sum(comp_product * H_full)
	total_vector = np.append(total_vector[:-1], H_c_all)
	
	normal_vec, d = compute_plane_4d(*realign_vectors)
	plane_energy = compute_energy_from_plane(normal_vec, d, total_vector[:-1])
	return total_vector[-1] - plane_energy