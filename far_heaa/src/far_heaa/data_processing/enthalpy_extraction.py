class EnthalpyExtractions:
    """
    A class to handle the extraction of binary enthalpy values.

    This class provides methods to extract enthalpy values or interaction parameters (e.g., omega)
    from a dictionary that maps binary element pairs to their associated enthalpy values.

    Methods:
        extract_binary_enthalpy(binary_dict, ele_pair):
            Extracts the enthalpy value for a specific element pair from the provided dictionary.
    """

    @staticmethod
    def extract_binary_enthalpy(binary_dict: dict[str, float], ele_pair: str) -> float:
        """
        Extract the binary enthalpy value from a given dictionary.

        Args:
            binary_dict (dict[str, float]): A dictionary containing binary element pairs as keys
                                            and their corresponding enthalpy values/omegas (float).
            ele_pair (str): A string representing the element pair for which to extract the enthalpy.
                            Example format: 'Fe-Ni'.

        Returns:
            float: The enthalpy value associated with the specified element pair.

        Raises:
            KeyError: If the element pair is not found in the binary dictionary.

        Example::

            binary_dict = {'Fe-Ni': -1.5, 'Ni-Cu': -0.8}
            ele_pair = 'Fe-Ni'

            # Extracts the enthalpy for the Fe-Ni pair
            enthalpy = EnthalpyExtractions.extract_binary_enthalpy(binary_dict, ele_pair)
            # enthalpy = -1.5
        """
        return binary_dict[ele_pair]
