from typing import Union
from itertools import combinations, combinations_with_replacement
import itertools
import numpy as np


def create_mol_grid(n, grid_size):
    mol_grid = np.round(np.linspace(0.01, 0.99, grid_size), 3)
    mol_combs = [list(i) for i in list(itertools.permutations(mol_grid, n - 1)) if sum(list(i)) <= 1]

    total_mol_combs = [i + [abs(1 - sum(i))] for i in mol_combs] + [[1/n]*n]
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
    check_mol_combs = [sum(i) for i in renew_tot]
    assert sum(check_mol_combs) == len(renew_tot)

    return renew_tot

def create_multinary(
        element_list: list[str],
        no_comb: Union[list[int], str]
) -> dict[int: list[str]]:
    """
	function that creates an alloy depending on no_comb value.
	:param element_list: The list of element symbols to create combinations for
	:param no_comb: A list of number to create combinations for. For ex. [2,3] means create binary and ternary. Can take
	value "all" for combinations from 2-7.
	:return: a list of n-nary lists containing str of the form 'E1-E2-E3...'
	"""
    if no_comb == "all":
        no_comb = list(range(2, 8))

    result = {}
    for idx, comb in enumerate(no_comb):
        temp_joint = list(combinations(element_list, comb))
        s = '-'  # To join the element
        final_joint = sorted(
            {s.join(sorted(pair)) for pair in temp_joint}
        )
        if final_joint:  # just a measure
            result[comb] = final_joint

    return result
