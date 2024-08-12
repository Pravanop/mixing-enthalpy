import json
import pandas as pd


class DataUtils:
    """
    Class for data loading and manipulation
    Attributes
    None
    Methods
    load_csv: Loads data from csv file
    load_json: Loads data from json file
    extract_ele_list: Extracts element list from provided data folder
    extract_binaryEnthalpy: Extracts pairwise mixing enthalpy values from provided data folder
    """

    @staticmethod
    def load_csv(
            input_data_path: str,
            lattice: str,
            source: str
    ) -> pd.DataFrame:
        """
        Loads mixing enthalpy data from csv file. Is currently deprecated as all input files should be provided in
        JSON format.
        Args:
            input_data_path: The folder path for the data folder. Recommended to create a Data/input data folder.
            lattice: lattice currently working on. Accepted values are 'bcc', 'fcc', 'hcp'
            source: a unique source identifier. The naming conventions for the input data can be read in the read me file for this folder.
        Returns:
            A pandas data frame consisting the binary mixing enthalpy keys as 'A-B' (alphabetically sorted) with
        its value in ev/atom.
        """
        file_path = f"{input_data_path}/{lattice}_{source}.csv"
        df = pd.read_csv(file_path)
        df.set_index('ele', inplace=True)
        return df

    @staticmethod
    def load_json(
            folder_path: str,
            lattice: str,
            source: str
    ) -> dict:
        """
        Loads binary mixing enthalpies from a json file for a particular source and a particular lattice.
        Args:
            folder_path: The folder path for the data folder. Recommended to create a Data/input data folder.
            lattice: lattice currently working on. Accepted values are 'bcc', 'fcc', 'hcp'
            source: a unique source identifier. The naming conventions for the input data can be read in the read me file for this folder.

        Returns:
            A dictionary consisting the binary mixing enthalpy keys as 'A-B' (alphabetically sorted) with
        its value in ev/atom.
        """
        with open(f'{folder_path}/{source}/{lattice}_{source}.json', 'r') as f:
            return json.load(f)

    @staticmethod
    def extract_ele_list(
            folder_path: str,
            lattice: str,
            source: str
    ) -> list[str]:
        """
        Extracts the element list from the provided text file in the input data folder.
        Note - The element text file should be comma delimited without any spaces between the commas and the element names. Can prove to be catastrophic if done so.
        Args:
           folder_path: The folder path for the data folder. Recommended to create a Data/input data folder.
            lattice: lattice currently working on. Accepted values are 'bcc', 'fcc', 'hcp'
            source: a unique source identifier. The naming conventions for the input data can be read in the read me file for this folder.

        Returns:
        A list of strings containing the element list from the provided text file.
        """
        with open(f"{folder_path}/{source}/element_list_{lattice}_{source}.txt", 'r') as f:
            ele_list = f.read()

        return ele_list.split(',')

    @staticmethod
    def dict_to_csv(
            results_dict: dict,
            source: str,
            lattice: str,
            n_alloy: int,
            output_folder: str,
    ) -> None:
        """
        Saves calculated mixing enthalpy values to a csv file in a specified folder. Is currently deprecated.
        Args:
            results_dict: The dictionary containing the calculated mixing enthalpy values. Should be of the form
                        {'A-B-C' : {
                        'mixing_enthalpy': <float>,
                        'config_entropy':<float>
                        },...}
                        The results_dict should ideally contains alloys of one particular dimensionality
            source: a unique source identifier. The naming conventions for the input data can be read in the read me file for this folder.
            lattice: lattice currently working on. Accepted values are 'bcc', 'fcc', 'hcp'
            n_alloy: the dimensionality of the alloys in the dict.
            output_folder: The folder path to store the csv file in

        Returns:
            None.
        """
        useful_names = {
            2: 'binary',
            3: 'tertiary',
            4: 'quaternary',
            5: 'quinary',
            6: 'senary',
            7: 'septenary',
            8: 'octanery'
        }

        temp_df = pd.DataFrame.from_dict(data=results_dict).transpose()
        path = f"{output_folder}/{useful_names[n_alloy]}_{lattice}_{source}.csv"
        temp_df.to_csv(path)

    @staticmethod
    def extract_binaryEnthalpy(
            binary_dict: dict,
            ele_pair: str
    ) -> float:
        """
        Finds the binary mixing enthalpy value from a dictionary for an unique key.
        Args:
            binary_dict: The dictionary of binary mixing enthalpies that are to be searched.
            ele_pair: The key in the form of 'A-B' (alphabetically sorted).

        Returns:
            Mixing enthalpy value in eV/atom.
            Raises KeyError if an element pair is not found.
        """
        return binary_dict[ele_pair]
