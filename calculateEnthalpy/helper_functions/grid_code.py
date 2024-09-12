from typing import Union
from itertools import combinations
import itertools
import numpy as np

"""
Two grid functions that allow for accessing 
The composition space - create_mol_grid, and
The phase space - create_multinary
"""


# def create_mol_grid(n: int,
# 					grid_size: int) -> np.array:
# 	"""
# 	Creates a grid of compositions between (0, 1) for n dimensions.
#
# 	Does this by finding all permutations of a given n for a numpy array ranging from 0 to 1. Then it filters
# 	through each permutation and checks if the sum is 1 as per the mol fraction constraint.
# 	Args:
# 		n: The dimensionality of the grid.
# 		grid_size: The number of points between 0 and 1 to be regarded. Values between 15-30 work best.
#
# 	Returns: A  numpy array of grid points of the form [[x1,y1,z1], [x2,y2,z2]...]
#
# 	"""
# 	mol_grid = np.round(np.linspace(0.01, 0.99, grid_size), 3)
# 	mol_combs = [list(i) for i in list(itertools.permutations(mol_grid, n - 1)) if sum(list(i)) <= 1]
#
# 	total_mol_combs = [i + [abs(1 - sum(i))] for i in mol_combs] + [[1 / n] * n]
# 	tot = []
# 	for i in total_mol_combs:
# 		tot.append(i)
# 		for j in list(itertools.permutations(i, len(i))):
# 			tot.append(list(j))
# 	counter_dict = {}
# 	renew_tot = []
# 	for i in tot:
# 		if str(i) in counter_dict:
# 			continue
# 		else:
# 			renew_tot.append(i)
# 			counter_dict[str(i)] = 1
#
# 	# ensure that the ratios chosen sum to one.
# 	check_mol_combs = [sum(i) for i in renew_tot]
# 	assert sum(check_mol_combs) == len(renew_tot) #since all are ones, the sum of sums of the grid should equal the total number of points.
#
# 	return np.array(renew_tot)

def create_mol_grid(n: int,
				grid_size: int):
	"""
	generate N-1 dimensional grid with a linear edge density of Npoints

	Parameters
	----------
	n : int > 1
		number of elements in alloy system.
	grid_size : int > 1
		linear point density on edges of grid.

	Returns
	-------
	out : list of vectors
		grid of compositions.

	"""
	N = n
	Npoints = grid_size
	start = 1e-8
	step = (1.0 - 2 * start) / (Npoints - 1)
	stop = 1 + step - start

	start = 0.0
	step = (1.0) / (Npoints - 1)
	stop = 1 + step

	x = np.arange(start, stop, step)
	xi = []
	for i in range(N - 1):
		xi.append(x)
	mesh = np.meshgrid(*xi)
	M = len(mesh)
	Xvecs = []
	for X in mesh:
		# print(X)
		Xvecs.append(np.reshape(X, X.size))

	points = []
	for i in range(len(Xvecs[0])):
		point = []
		for Xvec in Xvecs:
			point.append(Xvec[i])
		points.append(point)

	outpoints = []
	for point in points:
		isKeep = True
		cumsumx = 0
		for x in point:
			isKeep = isKeep and x <= 1 - cumsumx
			cumsumx += x  # -start/2.0
		if isKeep:
			outpoints.append(point)
	out = np.zeros((len(outpoints)+1, N))
	for i, pt in enumerate(outpoints):
		for j, xj in enumerate(pt):
			out[i, j] = xj
	out[:, N - 1] = 1 - np.sum(out[:, 0:(N - 1)], axis=1)
	out[-1] = np.array([1/N]*N)
	return np.round(out,4)
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
