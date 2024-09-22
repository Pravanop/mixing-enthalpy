from typing import Union, List, Dict, Literal
from itertools import combinations


class MultinaryCombinations:
    """
    A class to generate multinary combinations of elements.

    This class provides static methods to generate combinations of chemical elements
    based on the provided list of elements and desired combination sizes. The combinations
    can range from binary (2 elements) to septenary (7 elements), or be specified
    manually.

    Methods:
        create_multinary(element_list, no_comb, sort=True):
            Creates a dictionary of multinary combinations based on the desired
            combination sizes.

        generate_combinations(element_list, comb_size, sort=True):
            Generates specific combinations of elements for a given combination size.
    """

    @staticmethod
    def create_multinary(
        element_list: List[str],
        no_comb: Union[List[int], Literal["all"]],
        sort: bool = True,
    ) -> Dict[int, List[str]]:
        """
        Create multinary combinations of elements from the provided element list.

        Args:
                element_list (List[str]): A list of elements as strings (e.g., ['Fe', 'Ni', 'Co']).
                no_comb (Union[List[int], Literal['all']]): Either a list of integers representing the desired
                                                                                                        combination sizes (e.g., [2, 3, 4]) or the literal
                                                                                                        'all' to include combinations from binary (2) to septenary (7).
                sort (bool): Whether to sort the element combinations alphabetically. Default is True.

        Returns:
                Dict[int, List[str]]: A dictionary where the keys are the combination sizes and the values
                                                          are lists of element combinations (e.g., ['Co-Fe', 'Ni-Co']).
        """

        if no_comb == "all":
            no_comb = list(
                range(2, 8)
            )  # Default range from binary (2) to septenary (7)

        result = {}
        for comb_size in no_comb:
            result[comb_size] = MultinaryCombinations.generate_combinations(
                element_list, comb_size, sort
            )
        return result

    @staticmethod
    def generate_combinations(
        element_list: List[str], comb_size: int, sort: bool
    ) -> List[str]:
        """
        Generate element combinations of a specific size from the given element list.

        Args:
                element_list (List[str]): A list of elements as strings (e.g., ['Fe', 'Ni', 'Co']).
                comb_size (int): The size of the combinations to generate (e.g., 2 for binary combinations).
                sort (bool): Whether to sort each combination alphabetically.

        Returns:
                List[str]: A list of element combinations in string format (e.g., ['Co-Fe', 'Ni-Co']).
        """

        temp_combinations = list(combinations(element_list, comb_size))
        separator = "-"  # Standard separator for element combinations

        if sort:
            return sorted({separator.join(sorted(pair)) for pair in temp_combinations})

        return list({separator.join(pair) for pair in temp_combinations})
