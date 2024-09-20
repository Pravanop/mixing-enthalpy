class EnthalpyExtractions:
    @staticmethod
    def extract_binary_enthalpy(binary_dict: dict[str, float],
                                ele_pair: str) -> float:
        """
        Extract the binary enthalpy value from a given dictionary.

        Args:
            binary_dict (dict[str, float]): A dictionary containing binary pairs as keys
                                             and their corresponding enthalpy values/omegas as float.
            ele_pair (str): A string representing the element pair for which to extract the enthalpy.

        Returns:
            float: The enthalpy value associated with the specified element pair.

        Raises:
            KeyError: If the element pair is not found in the binary dictionary.
        """
        return binary_dict[ele_pair]

