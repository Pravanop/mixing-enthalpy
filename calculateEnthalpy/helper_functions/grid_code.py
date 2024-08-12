from typing import Union
from itertools import combinations
import itertools
import numpy as np

"""
Two grid functions that allow for accessing 
The composition space - create_mol_grid, and
The phase space - create_multinary
"""


def create_mol_grid(n: int,
					grid_size: int) -> np.array:
	"""
	Creates a grid of compositions between (0, 1) for n dimensions.

	Does this by finding all permutations of a given n for a numpy array ranging from 0 to 1. Then it filters
	through each permutation and checks if the sum is 1 as per the mol fraction constraint.
	Args:
		n: The dimensionality of the grid.
		grid_size: The number of points between 0 and 1 to be regarded. Values between 15-30 work best.

	Returns: A  numpy array of grid points of the form [[x1,y1,z1], [x2,y2,z2]...]

	"""
	mol_grid = np.round(np.linspace(0.01, 0.99, grid_size), 3)
	mol_combs = [list(i) for i in list(itertools.permutations(mol_grid, n - 1)) if sum(list(i)) <= 1]

	total_mol_combs = [i + [abs(1 - sum(i))] for i in mol_combs] + [[1 / n] * n]
	tot = []
	for i in total_mol_combs:
		tot.append(i)
		for j in list(itertools.permutations(i, len(i))):
			tot.append(list(j))
	counter_dict = {}
	renew_tot = []
	for i in tot:
		if str(i) in counter_dict:
			continue
		else:
			renew_tot.append(i)
			counter_dict[str(i)] = 1

	# ensure that the ratios chosen sum to one.
	check_mol_combs = [sum(i) for i in renew_tot]
	assert sum(check_mol_combs) == len(renew_tot) #since all are ones, the sum of sums of the grid should equal the total number of points.

	return np.array(renew_tot)


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
		)
		if final_joint:  # just a measure
			result[comb] = final_joint

	return result
