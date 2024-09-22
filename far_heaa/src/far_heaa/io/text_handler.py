from typing import List


class TextHandler:
    """
    A class to handle text file operations, such as extracting an element list from a text file.

    Methods:
        extract_ele_list(folder_path, lattice, source):
            Extracts the element list from a text file in the input data folder.
    """

    @staticmethod
    def extract_ele_list(folder_path: str, lattice: str, source: str) -> List[str]:
        """
        Extracts the element list from a text file in the specified folder.

        Args:
            folder_path (str): The path to the folder where the text file is located.
            lattice (str): The lattice type (e.g., "FCC", "BCC") for the element list.
            source (str): The source or name to distinguish between different datasets.

        Returns:
            List[str]: A list of elements extracted from the text file.

        Example:
            elements = TextHandler.extract_ele_list(folder_path='/data/', lattice='FCC', source='thermo_data')
            # This will read the file '/data/thermo_data/element_list_FCC_thermo_data.txt'
            # and return the elements as a list of strings.
        """
        with open(
            f"{folder_path}/{source}/element_list_{lattice}_{source}.txt", "r"
        ) as f:
            ele_list = f.read()

        return ele_list.split(",")
