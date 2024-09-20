class TextHandler:
    @staticmethod
    def extract_ele_list(folder_path: str, lattice: str, source: str) -> list[str]:
        """
        Extracts the element list from a text file in the input data folder.
        """
        with open(f"{folder_path}/{source}/element_list_{lattice}_{source}.txt", 'r') as f:
            ele_list = f.read()

        return ele_list.split(',')
