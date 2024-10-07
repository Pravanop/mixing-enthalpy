from pymatgen.core import Element
import numpy as np


def get_atomic_radii(elements: list[str]) -> np.array:
	"""
	Get atomic radii for a list of elements.

	Args:
		elements (list): List of elements.

	Returns:
		list: List of atomic radii.
	"""
	radii = []
	for element in elements:
		radii.append(Element(element).atomic_radius_calculated)
	return np.array(radii)


def get_electronegativity(elements: list[str]) -> np.array:
	"""
	Get electronegativity for a list of elements.

	Args:
		elements (list): List of elements.

	Returns:
		list: List of electronegativity.
	"""
	electronegativity = []
	for element in elements:
		electronegativity.append(Element(element).X)
	return np.array(electronegativity)


def get_ionization_energy(elements: list[str]) -> np.array:
	"""
	Get ionization energy for a list of elements.

	Args:
		elements (list): List of elements.

	Returns:
		list: List of ionization energy.
	"""
	ionization_energy = []
	for element in elements:
		ionization_energy.append(Element(element).ionization_energy)
	return np.array(ionization_energy)


def row_number(elements: list[str]) -> np.array:
	"""
	Get row number for a list of elements.

	Args:
		elements (list): List of elements.

	Returns:
		list: List of row numbers.
	"""
	row_number = []
	for element in elements:
		row_number.append(Element(element).row)
	return np.array(row_number)


def column_number(elements: list[str]) -> np.array:
	"""
	Get column number for a list of elements.

	Args:
		elements (list): List of elements.

	Returns:
		list: List of column numbers.
	"""
	column_number = []
	for element in elements:
		column_number.append(Element(element).group)
	return np.array(column_number)


def get_valence_electrons(elements: list[str]) -> np.array:
	"""
	Get valence electrons for a list of elements.

	Args:
		elements (list): List of elements.

	Returns:
		list: List of valence electrons.
	"""
	valence_electrons = []
	for element in elements:
		valence_electrons.append(Element(element).full_electronic_structure[-1][-1])
	return np.array(valence_electrons)


def get_atomic_number(elements: list[str]) -> np.array:
	"""
	Get atomic number for a list of elements.

	Args:
		elements (list): List of elements.

	Returns:
		list: List of atomic numbers.
	"""
	atomic_number = []
	for element in elements:
		atomic_number.append(Element(element).Z)
	return np.array(atomic_number)


def is_something(elements: list[str]) -> np.array:
	"""
	Check if element is a metal.

	Args:
		elements (list): List of elements.

	Returns:
		list: List of 1s and 0s.
	"""
	is_something_list = []
	for element in elements:
		ele = Element(element)
		temp = [int(ele.is_metalloid),
				int(ele.is_transition_metal),
				int(ele.is_post_transition_metal),
				int(ele.is_metal),
				int(ele.is_alkaline)]
		is_something_list.append(temp)
	return np.array(is_something_list)


def featurizer(elements: list[str]) -> np.array:
	features = [
		get_atomic_radii(elements),
		get_electronegativity(elements),
		get_ionization_energy(elements),
		row_number(elements),
		column_number(elements),
		get_valence_electrons(elements),
		get_atomic_number(elements)]
	
	features_array = np.zeros((len(elements), len(features)))
	for i in range(len(features)):
		features_array[:, i] = features[i]
	
	is_something_feature = is_something(elements)
	features_array = np.concatenate((features_array, is_something_feature), axis=1)
	features_desc = ['Atomic Radii', 'Electronegativity', 'Ionization Energy', 'Row Number', 'Column Number',
					 'Valence Electrons', 'Atomic Number', 'Is Metalloid', 'Is Transition',
					 'Is Post Transition', 'Is Metal', 'Is Alkaline']
	return features_array, features_desc
