from typing import Union
from itertools import combinations
def create_multinary(
		element_list: list[str] ,
		no_comb: Union[list[int] , str]
		) -> list[list[str]] :
	"""
	function that creates an alloy depending on no_comb value.
	:param element_list: The list of element symbols to create combinations for
	:param no_comb: A list of number to create combinations for. For ex. [2,3] means create binary, ternary. Can take
	value "all" for combinations from 2-7.
	:return: a list of n-nary lists containing str of the form 'E1-E2-E3...'
	"""
	if no_comb == "all" :
		no_comb = list(range(2 , 8))
	
	result = []
	for idx , comb in enumerate(no_comb) :
		temp_joint = list(combinations(element_list , comb))
		s = '-'  # To join the element
		final_joint = sorted(
				{s.join(sorted(pair)) for pair in temp_joint}
				)
		if final_joint :  # just a measure
			result.append(final_joint)
	
	return result

def create_multinary_unsorted(
		element_list: list[str] ,
		no_comb: Union[list[int] , str]
		) -> list[list[str]] :
	"""
	function that creates an alloy depending on no_comb value.
	:param element_list: The list of element symbols to create combinations for
	:param no_comb: A list of number to create combinations for. For ex. [2,3] means create binary, ternary. Can take
	value "all" for combinations from 2-7.
	:return: a list of n-nary lists containing str of the form 'E1-E2-E3...'
	"""
	if no_comb == "all" :
		no_comb = list(range(2 , 8))
	
	result = []
	for idx , comb in enumerate(no_comb) :
		temp_joint = list(combinations(element_list , comb))
		s = '-'  # To join the element
		final_joint = sorted(
				{s.join(pair) for pair in temp_joint}
				)
		if final_joint :  # just a measure
			result.append(final_joint)
	
	return result