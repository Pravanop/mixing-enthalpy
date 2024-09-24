import os.path
from typing import List


class DirHandler:
    """
    A class to handle directory creation, including creating nested directories recursively.

    Methods:
            mkdir_recursive(folders, folder_path):
                    Creates nested directories in the specified path if they do not already exist.
    """

    @staticmethod
    def mkdir_recursive(folders: List[str], folder_path: str) -> str:
        """
        Recursively creates nested directories in the specified folder path.

        Args:
                folders (List[str]): A list of folder names to create in the specified path.
                folder_path (str): The base path where the folders should be created.

        Returns:
                str: The final path after all directories have been created.

        Example::
        
                DirHandler.mkdir_recursive(folders=['data', 'results', 'plots'], folder_path='/home/user/project/')
                # This will create '/home/user/project/data/results/plots/' if it does not exist.
        """

        # Ensure the folder path ends with a '/'
        if "/" != folder_path[-1]:
            folder_path += "/"
        if not os.path.exists(folder_path):
            os.mkdir(folder_path)

        # Iterate through each folder in the list and create the path
        for folder in folders:
            folder_path = folder_path + folder + "/"
            if not os.path.exists(folder_path):
                os.mkdir(folder_path)

        return folder_path
