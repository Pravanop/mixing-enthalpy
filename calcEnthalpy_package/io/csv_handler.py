import pandas as pd

class CSVHandler:
    @staticmethod
    def load_csv(input_data_path: str, lattice: str, source: str) -> pd.DataFrame:
        """
        Loads mixing enthalpy data from csv file.
        """
        file_path = f"{input_data_path}/{lattice}_{source}.csv"
        df = pd.read_csv(file_path)
        df.set_index('ele', inplace=True)
        return df

    @staticmethod
    def dict_to_csv(results_dict: dict, source: str, lattice: str, n_alloy: int, output_folder: str) -> None:
        """
        Saves calculated mixing enthalpy values to a csv file in a specified folder.
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
