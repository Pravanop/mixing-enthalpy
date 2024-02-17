"""
Gets binary pairs for SQS creation
"""
from calculateEnthalpy.create_alloy_comp import create_multinary_unsorted

def get_binary_pairs(ele_list: list) -> list[str]:
	"""
	
	:param ele_list:
	:return:
	"""
	return create_multinary_unsorted(
		element_list = ele_list ,
		no_comb = [2]
		)[0]
