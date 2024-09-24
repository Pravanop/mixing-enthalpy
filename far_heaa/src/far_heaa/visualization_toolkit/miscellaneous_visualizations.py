import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pymatgen.analysis.phase_diagram import PDPlotter
import matplotlib as mpl
from typing import List, Tuple, Literal

from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.io.dir_handler import DirHandler
from far_heaa.io.json_handler import JSONHandler
from far_heaa.math_operations.thermo_calculations import ThermoMaths
from far_heaa.phase_diagram.grid_iterators import GridIterator


class MiscellaneousVisualizations:
    """
    A class for generating various visualizations related to alloy compositions, including convex hull plots
    and heatmaps for binary mixing enthalpies.

    Args:
            meta_data (dict): A dictionary containing metadata such as grid size, file paths, and flags.
            save_flag (bool): A flag indicating whether to save the visualizations.

    Attributes:
            grid_iterator (GridIterator): An instance of the GridIterator for iterating through grid data.
            tm (ThermoMaths): An instance of the ThermoMaths class for thermodynamic calculations.
            save_flag (bool): A flag indicating whether to save the generated plots.
    Methods:
            plot_convex_hull: Plots the convex hull for a given alloy composition at a specified temperature.
            heatmap: Creates a heatmap of binary mixing enthalpy for a list of elements.
    """

    def __init__(self, meta_data: dict, save_flag: bool):
        """
        Initializes the MiscellaneousVisualizations class with the given lattice structure and metadata.

        Args:
                meta_data (dict): A dictionary containing metadata such as grid size, file paths, and flags.
                save_flag (bool): A flag indicating whether to save the visualizations.
        """

        self.mol_grid_size = 40
        self.flags = meta_data["flags"]
        self.tm = ThermoMaths()
        grid_size = meta_data["grid_size"]
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
            grid_size=grid_size,
            tm=self.tm,
            data=data,
            end_member=end_member,
            api_key=meta_data["api_key"],
            flags=meta_data["flags"],
        )
        self.save_flag = save_flag

    def plot_convex_hull(
        self, composition: List[str], temperature: float
    ) -> Tuple[plt.Axes, plt.Figure]:
        """
        Plots the convex hull for a given alloy composition at a specified temperature.

        Args:
                composition (List[str]): A list of elements in the alloy composition.
                temperature (float): The temperature in Kelvin.

        Returns:
                Tuple[plt.Axes, plt.Figure]: A tuple containing the matplotlib axes and figure objects for the plot.

        Example:
                ax, fig = self.plot_convex_hull(composition=['Fe', 'Ni'], temperature=1200)
        """
        conv_hull = self.grid_iterator.convex_hull.make_convex_hull(
            temperature=temperature, composition=composition, batch_tag=False
        )
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        mpl.rcParams.update({"font.size": 10})

        PDPlotter(
            phasediagram=conv_hull, show_unstable=0.02, backend="matplotlib"
        ).get_plot(
            label_unstable=False,
            label_stable=False,
            process_attributes=True,
            ax=ax,
            fill=False,
        )

        if len(composition) == 2:
            ax.set_xlim([0, 1])
            ax.set_xlabel("x", fontsize=10)
            ax.set_ylabel("E hull (eV/atom)", fontsize=10)
            plt.subplots_adjust(bottom=0.12, left=0.15, right=0.95, top=0.85)
            ax.set_title(f'{"".join(sorted(composition, reverse=True))}')

        if len(composition) == 3:
            plt.subplots_adjust(top=0.95)
            vert = [(0.5, 0.9), (-0.05, 0), (1.02, 0)]
            composition.reverse()
            for idx, i in enumerate(composition):
                ax.text(s=i, x=vert[idx][0], y=vert[idx][1], fontsize=12)

        if self.save_flag:
            updated_folder_path = DirHandler.mkdir_recursive(
                folders=[
                    "convex_hull",
                    f"{len(composition)}",
                    f'{"".join(sorted(composition))}',
                ],
                folder_path="../plots",
            )
            save_path = (
                f"{updated_folder_path}/convex_hull_{temperature}_equi.png"
                if self.flags["equi_flag"]
                else f"{updated_folder_path}/convex_hull_{temperature}.png"
            )
            fig.savefig(fname=save_path, dpi=100)

        return ax, fig

    def heatmap(
        self, element_list: List[str], lattice: Literal["FCC", "BCC", "HCP"], sns=None
    ) -> plt.Figure:
        """
        Creates a heatmap of binary mixing enthalpy for a list of elements.

        Args:
                element_list (List[str]): A list of elements for which to create the heatmap.
                lattice (Literal['FCC', 'BCC', 'HCP']): The lattice to display (e.g., 'FCC', 'BCC').
                sns (optional): Seaborn instance for customized plotting.

        Returns:
                plt.Figure: The matplotlib figure object for the heatmap.

        Example:
                heatmap = self.heatmap(element_list=['Fe', 'Ni', 'Al'], genre='unbiased')
        """
        binaries = list(
            MultinaryCombinations.create_multinary(element_list, no_comb=[2])[2]
        )
        final_data = np.empty((len(element_list), len(element_list)))
        final_data[:] = np.nan

        for pair in binaries:
            pair_list = pair.split("-")
            value = (
                int(self.grid_iterator.data[pair][lattice] * 1000 / 4)
                if self.flags["correction"]
                else int(self.grid_iterator.data[pair][lattice] * 1000)
            )
            final_data[element_list.index(pair_list[0])][
                element_list.index(pair_list[1])
            ] = value
            final_data[element_list.index(pair_list[1])][
                element_list.index(pair_list[0])
            ] = value

        alloy = pd.DataFrame(final_data, index=element_list, columns=element_list)
        fig = sns.heatmap(
            alloy,
            cmap="coolwarm",
            annot=True,
            fmt="g",
            mask=alloy.isnull(),
            square=True,
            cbar_kws={"label": "Binary H$_{mix}$ (meV/atom)"},
            linewidths=2,
            linecolor="white",
        )

        return fig
