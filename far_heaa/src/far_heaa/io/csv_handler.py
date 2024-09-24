import pandas as pd
from typing import Dict


class CSVHandler:
    """
    A class for handling CSV file operations such as loading data from and saving data to CSV files.

    Methods:
        load_csv(input_data_path, lattice, source):
            Loads mixing enthalpy data from a CSV file.

        dict_to_csv(results_dict, source, lattice, n_alloy, output_folder):
            Saves calculated mixing enthalpy values from a dictionary to a CSV file.
    """

    @staticmethod
    def load_csv(input_data_path: str, lattice: str, source: str) -> pd.DataFrame:
        """
        Loads mixing enthalpy data from a CSV file.

        Args:
            input_data_path (str): The path to the folder containing the CSV file.
            lattice (str): The lattice type (e.g., "FCC", "BCC") to load the corresponding CSV file.
            source (str): The source of the data to distinguish between different datasets.

        Returns:
            pd.DataFrame: A DataFrame containing the loaded mixing enthalpy data with elements as the index.

        Example::
        
            df = CSVHandler.load_csv(input_data_path='/data/', lattice='FCC', source='thermo_data')
            # Loads the data from '/data/FCC_thermo_data.csv' into a DataFrame.
        """
        file_path = f"{input_data_path}/{lattice}_{source}.csv"
        df = pd.read_csv(file_path)
        df.set_index("ele", inplace=True)
        return df

    @staticmethod
    def dict_to_csv(
        results_dict: Dict, source: str, lattice: str, n_alloy: int, output_folder: str
    ) -> None:
        """
        Saves calculated mixing enthalpy values to a CSV file in a specified folder.

        Args:
            results_dict (Dict): A dictionary containing calculated mixing enthalpy values to be saved.
            source (str): The source or name of the data being saved (e.g., "thermo_data").
            lattice (str): The lattice type (e.g., "FCC", "BCC") used in the calculations.
            n_alloy (int): The number of elements in the alloy (e.g., 2 for binary, 3 for tertiary).
            output_folder (str): The path to the folder where the CSV file will be saved.

        Returns:
            None

        Example::
        
            CSVHandler.dict_to_csv(results_dict, source='thermo_data', lattice='FCC', n_alloy=3, output_folder='/results/')
            # Saves the data to '/results/tertiary_FCC_thermo_data.csv'.
        """
        useful_names = {
            2: "binary",
            3: "tertiary",
            4: "quaternary",
            5: "quinary",
            6: "senary",
            7: "septenary",
            8: "octanery",
        }
        temp_df = pd.DataFrame.from_dict(data=results_dict).transpose()
        path = f"{output_folder}/{useful_names[n_alloy]}_{lattice}_{source}.csv"
        temp_df.to_csv(path)
