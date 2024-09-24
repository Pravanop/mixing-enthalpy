from far_heaa.io.json_handler import JSONHandler


class MetadataHandler:
    """
    A class to handle metadata operations, including loading, accessing, and updating metadata stored in a JSON file.
    It is recommended to use this class for everything. One can get all the keys from access_metadata_keys method
    Args:
            folder_path (str, optional): The folder path where the metadata JSON file is located. Defaults to '../database/'.
            file_name (str, optional): The name of the metadata JSON file (without extension). Defaults to 'metadata'.

    Attributes:
            meta_data (dict): The loaded metadata from the JSON file.
    """

    def __init__(self, folder_path: str = "../database/", file_name: str = "metadata"):
        """
        Initializes the MetadataHandler class by loading the metadata JSON file.

        Args:
                folder_path (str, optional): The folder path where the metadata JSON file is located. Defaults to '../database/'.
                file_name (str, optional): The name of the metadata JSON file (without extension). Defaults to 'metadata'.
        """
        self.meta_data = JSONHandler.load_json(
            file_name=file_name, folder_path=folder_path
        )
        for flag, value in self.meta_data["flags"].items():
            if value == "True":
                self.meta_data["flags"][flag] = True
            elif value == "False":
                self.meta_data["flags"][flag] = False

    @property
    def get_metadata(self) -> dict:
        """
        Retrieves the metadata as a dictionary.

        Returns:
                dict: The loaded metadata.
        """
        return self.meta_data

    def access_metadata_keys(self) -> dict.keys:
        """
        Accesses the keys of the metadata dictionary.

        Returns:
                dict_keys: The keys of the metadata.
        """

        self._print_nested_keys(self.meta_data)

    def _print_nested_keys(self, d: dict, indent: int = 0) -> None:
        """
        Recursively prints all the keys in a nested dictionary with indentation for each level.

        Args:
                d (dict): The dictionary whose keys need to be printed.
                indent (int, optional): The indentation level for nested keys. Defaults to 0.

        Example::
        
                print_nested_keys(nested_dict)
        """
        for key, value in d.items():
            print("    " * indent + str(key))  # Print the key with indentation
            if isinstance(
                value, dict
            ):  # If the value is a dictionary, recursively print its keys
                self._print_nested_keys(value, indent + 1)

    def update_metadata(self, key: str, value) -> None:
        """
        Updates the value of a specific key in the metadata.

        Args:
                key (str): The key in the metadata dictionary to update.
                value: The new value to assign to the key.

        Raises:
                KeyError: If the provided key is not present in the metadata.
        """
        if key in self.meta_data:
            self.meta_data[key] = value
        else:
            raise KeyError(f"Metadata key {key} not present")
