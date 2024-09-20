from typing import Union, List, Dict
from itertools import combinations


class MultinaryCombinations:
	"""
	A class to generate combinations of elements for alloy systems.

	Methods
	-------
	create_multinary(element_list: List[str], no_comb: Union[List[int], str], sort: bool = True) -> Dict[int, List[str]]:
		Generates combinations of elements for specified combination sizes.
	generate_combinations(element_list: List[str], comb_size: int, sort: bool) -> List[str]:
		Generates and optionally sorts combinations of elements.
	"""
	
	@staticmethod
	def create_multinary(
			element_list: List[str],
			no_comb: Union[List[int], str],
			sort: bool = True
	) -> Dict[int, List[str]]:
		"""
		Generates combinations of elements for specified combination sizes.

		Parameters
		----------
		element_list : list of str
			List of element symbols to create combinations for.
		no_comb : list of int or "all"
			Combination sizes to create. If "all", generates combinations from 2 to 7.
		sort : bool, optional
			Whether to sort elements within each combination (default is True).

		Returns
		-------
		dict of int : list of str
			Dictionary where keys are combination sizes and values are lists of combinations.
		"""
		if no_comb == "all":
			no_comb = list(range(2, 8))  # Default range from binary (2) to septenary (7)
		
		result = {}
		for comb_size in no_comb:
			result[comb_size] = MultinaryCombinations.generate_combinations(
				element_list, comb_size, sort
			)
		return result
	
	@staticmethod
	def generate_combinations(
			element_list: List[str],
			comb_size: int,
			sort: bool
	) -> List[str]:
		"""
		Generates and optionally sorts combinations of elements.

		Parameters
		----------
		element_list : list of str
			List of element symbols.
		comb_size : int
			The size of combinations to generate (e.g., binary, ternary).
		sort : bool
			Whether to sort elements within each combination (default is True).

		Returns
		-------
		list of str
			A list of generated combinations of elements.
		"""
		temp_combinations = list(combinations(element_list, comb_size))
		separator = '-'  # Standard separator for element combinations
		
		if sort:
			return sorted({separator.join(sorted(pair)) for pair in temp_combinations})
		else:
			return list({separator.join(pair) for pair in temp_combinations})
