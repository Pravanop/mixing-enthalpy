import json

class JSONHandler:
    @staticmethod
    def load_json(folder_path: str, file_name: str) -> dict:
        """
        Loads binary mixing enthalpies from a json file for a particular source and a particular lattice.
        """
        with open(f'{folder_path}/{file_name}.json', 'r') as f:
            return json.load(f)
