import json
from typing import Dict


class JSONHandler:
    """
    A class to handle reading and writing JSON files.

    Methods:
        load_json(folder_path, file_name):
            Loads data from a JSON file.
    """

    @staticmethod
    def load_json(folder_path: str, file_name: str) -> Dict:
        """
        Loads binary mixing enthalpies from a JSON file for a particular source and lattice.

        Args:
            folder_path (str): The path to the folder where the JSON file is located.
            file_name (str): The name of the JSON file (without the '.json' extension) to load.

        Returns:
            Dict: The contents of the JSON file as a dictionary.

        Example::
        
            data = JSONHandler.load_json(folder_path='/data/', file_name='mixing_enthalpy')
            # This will load the data from '/data/mixing_enthalpy.json' as a dictionary.
        """
        with open(f"{folder_path}/{file_name}.json", "r") as f:
            return json.load(f)
