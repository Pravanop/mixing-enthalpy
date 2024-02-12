from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure
def getAPIKey(api_key_file: str='callMpAPI/api_key.txt') -> str:
    """
    Gets API key for MP from a text file to be filled by user.
    :param api_key_file: str to file path, usually in same directory
    :return: str of APIkey
    """
    with open(api_key_file, 'r') as f:
        return f.read()

def convert_to_conventional(structure: Structure) -> Structure:
    """
    Given a primitive structure it converts to conventional setting
    :param structure: Pymatgen primitive structure
    :return: Pymatgen conventional structure
    """
    SGA = SpacegroupAnalyzer(structure)
    return SGA.get_conventional_standard_structure()