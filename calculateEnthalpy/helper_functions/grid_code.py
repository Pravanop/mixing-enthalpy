from typing import Union
from itertools import combinations
import numpy as np


def create_mol_grid(n: int, grid_size: int):
	"""
	Generate N-1 dimensional grid with a linear edge density of grid_size.

	Parameters
	----------
	n : int > 1
		Number of elements in alloy system.
	grid_size : int > 1
		Linear point density on edges of grid.

	Returns
	-------
	out : ndarray
		Grid of compositions (N points).
	"""
	if n < 2:
		raise ValueError("n must be greater than 1")
	
	# Create a grid in N-1 dimensions
	linspace = np.linspace(0, 1, grid_size)
	mesh = np.array(np.meshgrid(*[linspace] * (n - 1)))
	
	# Reshape the meshgrid to generate all possible combinations of mol fractions
	points = mesh.T.reshape(-1, n - 1)
	
	# Filter valid points where the sum of mol fractions <= 1
	valid_points = points[np.sum(points, axis=1) <= 1]
	
	# Compute the last component to ensure sum to 1
	last_column = 1 - np.sum(valid_points, axis=1, keepdims=True)
	
	# Combine valid points and last column
	out = np.hstack([valid_points, last_column])
	
	# Add a final point where all elements have equal mol fraction (homogeneous)
	homogeneous_point = np.full((1, n), 1 / n)
	out = np.vstack([out, homogeneous_point])
	
	return np.round(out, 4)


# return out

def create_multinary(
		element_list: list[str],
		no_comb: Union[list[int], str]
) -> dict[int: list[str]]:
	"""

	Args:
		element_list: The list of element symbols to create combinations for
		no_comb:  A list of number to create combinations for. For ex. [2,3] means create binary and ternary. Can take
	value "all" for combinations from 2-7.

	Returns: a dictionary of the form
	{
	3: ['E1-E2-E3', 'E1-E3-E4'....],
	2: ....
	}

	"""
	if no_comb == "all": #taking septenary as the limit, rarely used option.
		no_comb = list(range(2, 8))

	result = {}
	for idx, comb in enumerate(no_comb):
		temp_joint = list(combinations(element_list, comb))
		s = '-'  # To join the element, standard protocol
		final_joint = sorted(
			{s.join(sorted(pair)) for pair in temp_joint}
			# {s.join(pair) for pair in temp_joint}
		)
		if final_joint:  # just a measure
			result[comb] = final_joint

	return result
