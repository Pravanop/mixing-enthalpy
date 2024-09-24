from typing import Literal, List
from far_heaa.io.dir_handler import DirHandler
from far_heaa.math_operations.thermo_calculations import ThermoMaths
from far_heaa.phase_diagram.grid_iterators import GridIterator
from far_heaa.io.json_handler import JSONHandler
import matplotlib.pyplot as plt


class Visualizations:
    """
    A class used to represent visualizations for HEA.

    Attributes:
        lattice (Literal['FCC', 'BCC', 'HCP', 'min']): The lattice structure.
        mol_grid_size (int): The size of the mole fraction grid.
        tm (ThermoMaths): An instance of ThermoMaths for thermodynamic calculations.
        grid_size (int): The size of the grid.
        flags (dict): A dictionary of flags for various options.
        folder_path (str): The path to the folder where plots will be saved.
        grid_iterator (GridIterator): An instance of GridIterator for iterating over the grid.
     Methods:
         save_figure(folders: List[str], file_name: str, fig: plt.Figure) -> None:
         Saves the figure in the specified folders with different file name suffixes based on flags.
    """

    def __init__(self, lattice: Literal["FCC", "BCC", "HCP", "min"], meta_data: dict):
        """
        Initialize the Vizualization class with the lattice type, and metadata.

        Args:
                lattice (Literal['FCC', 'BCC', 'HCP', 'min']): The lattice structure. 'min' is whichever lattice gives minimum mixing enthalpy for that mole fraction.
                meta_data (dict): A dictionary containing metadata like grid size, file paths, and flags.

        """

        self.lattice = lattice
        self.mol_grid_size = 40
        self.tm = ThermoMaths()
        self.grid_size = meta_data["grid_size"]
        self.flags = meta_data["flags"]
        if meta_data["flags"]["correction"]:
            data = JSONHandler.load_json(
                folder_path=meta_data["folder_path"],
                file_name=meta_data["file_name"]["biased"],
            )
        else:
            data = JSONHandler.load_json(
                folder_path=meta_data["folder_path"],
                file_name=meta_data["file_name"]["unbiased"],
            )

        end_member = JSONHandler.load_json(
            folder_path=meta_data["folder_path"], file_name=meta_data["end_member"]
        )

        self.grid_iterator = GridIterator(
            grid_size=self.grid_size,
            tm=self.tm,
            data=data,
            end_member=end_member,
            api_key=meta_data["api_key"],
            flags=meta_data["flags"],
        )

        self.folder_path = "../plots"

    def save_figure(self, folders: List[str], file_name: str, fig: plt.Figure):
        """
        Saves the figure in the specified folders with different file name suffixes based on flags.

        Args:
            folders (list[str]): List of folder names to create recursively.
            file_name (str): The base name of the file to save the figure as.
            fig (plt.Figure): The matplotlib figure object to save.

        Returns:
            None
        """
        updated_folder_path = DirHandler.mkdir_recursive(
            folders=folders, folder_path=self.folder_path
        )

        if self.flags["equi_flag"] == True:
            fig.savefig(fname=f"{updated_folder_path}{file_name}_equi.png", dpi=100)

        if not self.flags["im_flag"]:
            fig.savefig(fname=f"{updated_folder_path}{file_name}_wo_im.png", dpi=100)

        if not self.flags["correction"]:
            fig.savefig(
                fname=f"{updated_folder_path}{file_name}_wo_correction.png", dpi=100
            )

        else:
            fig.savefig(fname=f"{updated_folder_path}{file_name}.png", dpi=100)
