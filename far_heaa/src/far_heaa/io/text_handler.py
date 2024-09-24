from typing import List


class TextHandler:
    """
    A class to handle text file operations, such as extracting an element list from a text file.

    Methods:
        extract_ele_list(folder_path, lattice, source):
            Extracts the element list from a text file in the input data folder.
    """

    @staticmethod
    def extract_ele_list(folder_path: str, file_name: str) -> List[str]:
        """
        Extracts the element list from a text file in the specified folder.

        Args:
            folder_path (str): The path to the folder where the text file is located.
            file_name (str): file name of the text file.

        Returns:
            List[str]: A list of elements extracted from the text file.

        """
        with open(f"{folder_path}/{file_name}.txt", "r") as f:
            ele_list = f.read()

        return ele_list.split(",")
