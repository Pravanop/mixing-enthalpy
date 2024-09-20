class EnthalpyExtractions:
    @staticmethod
    def extract_binary_enthalpy(binary_dict: dict, ele_pair: str) -> float:
        """
        Finds the binary mixing enthalpy value from a dictionary for a unique key.
        """
        return binary_dict[ele_pair]
