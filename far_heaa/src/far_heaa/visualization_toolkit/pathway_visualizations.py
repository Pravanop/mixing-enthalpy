import itertools
from typing import List, Dict, Tuple, Literal

import numpy as np
from matplotlib import pyplot as plt

from far_heaa.grids_and_combinations import grid_creation
from far_heaa.math_operations.thermo_calculations import ThermoMaths
from far_heaa.visualization_toolkit.visualizations import Visualizations


class PathwayVisualizations(Visualizations):
    """
    A class for visualizing reaction pathways for a given alloy composition and lattice structure.

    Args:
            composition (list): A list of elements in the alloy system.
            lattice (str): The lattice structure to consider (e.g., 'FCC', 'BCC').
            meta_data (dict): A dictionary containing metadata such as grid size, file paths, and flags.
            save_flag (bool): A flag indicating whether to save the visualizations.

    Attributes:
            tm (ThermoMaths): An instance of the ThermoMaths class for thermodynamic calculations.
            grid_iterator (GridIterator): An instance of the GridIterator class for grid-based iteration.
            temp_grid (list): A list of temperature values for the convex hull calculations.
            conv_hull (dict): A dictionary containing convex hull data for the alloy system.
            mol_grid_size (int): The size of the mole fraction grid.
            save_flag (bool): A flag indicating whether to save the generated plots.

    Methods:
            get_rP: Calculates reaction pathways for the alloy system and returns miscibility temperatures for each pathway.
            get_n_colors_from_cmap: Generates N evenly spaced colors from a given colormap.
            text_segregators: Groups text labels by their proximity and returns their adjusted positions for plotting.
            find_indices_in_range: Finds groups of indices in the input list where values are within the specified threshold.
            plot_rP: Plots the reaction pathways for the alloy system and the corresponding miscibility temperatures.

    """

    def __init__(
        self,
        composition: List[str],
        lattice: Literal["FCC", "BCC", "HCP", "min"],
        meta_data: dict,
        save_flag: bool,
    ):
        """
        Initializes the PathwayVisualizations class with the alloy composition, lattice type, and metadata.

        Args:
                composition (list): A list of elements in the alloy system.
                lattice (Literal['FCC', 'BCC', 'HCP', 'min']): The lattice structure. 'min' is whichever lattice gives minimum mixing enthalpy for that mole fraction.
                meta_data (dict): A dictionary containing metadata such as grid size, file paths, and flags.
                save_flag (bool): A flag indicating whether to save the visualizations.
        """

        super().__init__(lattice=lattice, meta_data=meta_data)

        self.x = None
        self.composition = composition
        t_max = max([self.tm.avg_T_melt(i, mol_ratio=[]) for i in self.composition])
        self.temp_grid = list(np.linspace(0, t_max, 30))
        self.conv_hull = self.grid_iterator.temp_iterator(
            composition=self.composition, temp_grid=self.temp_grid
        )

        self.mol_grid_size = 5
        self.save_flag = save_flag

    def get_rP(self) -> Dict[str, List[np.ndarray]]:
        """
        Calculates reaction pathways for the alloy system and returns miscibility temperatures for each pathway.

        Returns:
                dict: A dictionary where the keys are reaction pathways (e.g., 'Fe-Ni') and the values are lists of
                miscibility temperatures along each pathway.

        Example::
        
                path_dict = self.get_rP()
        """
        all_pathways = list(
            itertools.permutations(self.composition, len(self.composition))
        )

        path_dict = {}
        self.x = list(np.linspace(0, 1, self.mol_grid_size))

        for path in all_pathways:
            path = list(path)
            count = 0
            temp_path = []

            while count + 1 < len(path):
                start = path[: count + 1]
                end = [path[count + 1]]
                total = start + end
                mol_grid = grid_creation.CompositionGrid.create_high_sym_mol_grid(
                    change_idx=list(range(len(total) - 1)),
                    x=self.x,
                    n=len(total),
                    N=len(start),
                )
                mol_grid = mol_grid[::-1]
                _, misc_temp = self.grid_iterator.misc_temperature_across_grid(
                    composition=total,
                    mol_grid_size=mol_grid,
                    lattice=self.lattice,
                    phase_flag=True,
                    conv_hull=self.conv_hull,
                    temp_grid=self.temp_grid,
                )
                for idx, temp in enumerate(misc_temp):
                    if temp == -1:
                        misc_temp[idx] = 5000

                if count == 0:
                    misc_temp[0] = self.tm.avg_T_melt(
                        composition=path[: count + 1][0], mol_ratio=[]
                    )

                temp_path.append(misc_temp)
                count += 1
            path_dict["-".join(path)] = temp_path

        return path_dict

    @staticmethod
    def get_n_colors_from_cmap(cmap_name: str, N: int) -> List[tuple]:
        """
        Generates N evenly spaced colors from a given colormap.

        Args:
                cmap_name (str): The name of the colormap to use.
                N (int): The number of colors to generate.

        Returns:
                List[tuple]: A list of RGB tuples representing the colors.
        """
        cmap = plt.get_cmap(cmap_name)
        colors = [cmap(i) for i in np.linspace(0, 1, N)]
        return colors

    def text_segregators(self, texts: np.ndarray) -> dict[str, list[float]]:
        """
        Groups text labels by their proximity and returns their adjusted positions for plotting.

        Args:
                texts (np.ndarray): An array of text labels and their positions.

        Returns:
                dict: A dictionary where the keys are concatenated text labels and the values are their adjusted positions.
        """
        x_dict = {}
        for text in texts:
            if text[1] not in x_dict:
                x_dict[text[1]] = [text]
            else:
                x_dict[text[1]].append(text)

        text_dict = {}
        for key, value in x_dict.items():
            value = np.array(value)
            y = value[:, -1].astype(float)
            list_together = self.find_indices_in_range(y, threshold=200)
            for i in list_together:
                temp_text = "".join(f"{value[j, 0]}, " for j in i)[:-2]
                text_dict[temp_text] = [
                    float(value[i[0], 1]) - 0.02 * len(temp_text),
                    float(value[i[0], 2]),
                ]

        return text_dict

    @staticmethod
    def find_indices_in_range(
        float_list: np.ndarray, threshold: float = 200
    ) -> List[List[int]]:
        """
        Finds groups of indices in the input list where values are within the specified threshold.

        Args:
                float_list (np.ndarray): A list of float values to compare.
                threshold (float, optional): The maximum difference between values to group them. Default is 200.

        Returns:
                List[List[int]]: A list of lists where each inner list contains indices of grouped values.
        """
        indices_groups = []
        visited = set()

        for i in range(len(float_list)):
            if i in visited:
                continue

            current_group = [i]
            for j in range(i + 1, len(float_list)):
                if abs(float_list[i] - float_list[j]) <= threshold:
                    current_group.append(j)
                    visited.add(j)

            indices_groups.append(current_group)

        return indices_groups

    def plot_rP(self) -> Tuple[plt.Axes, plt.Figure]:
        """
        Plots the reaction pathways for the alloy system and the corresponding miscibility temperatures.

        Returns:
                Tuple[plt.Axes, plt.Figure]: A tuple containing the matplotlib axes and figure objects for the plot.

        Example::
        
                ax, fig = self.plot_rP()
        """
        path_dict = self.get_rP()
        cmap = self.get_n_colors_from_cmap(
            cmap_name="coolwarm", N=len(self.composition)
        )
        cmap_dict = dict(zip(self.composition, cmap))

        texts = []
        fig, ax = plt.subplots()

        for key, value in path_dict.items():
            path = key.split("-")
            for idx, i in enumerate(value):
                x = np.array(self.x) + idx
                y = np.array(i)
                ax.axvline(x=x[0], color="#BBBBBB", alpha=0.3, linestyle="--", zorder=0)
                temp = key.split("-")[: idx + 1]
                texts.append(["".join(sorted(temp)), x[0], y[0]])

                if idx == len(value) - 1:
                    texts.append(["".join(sorted(self.composition)), x[-1], y[-1]])

                ax.plot(x, y, zorder=0, color=cmap_dict[path[0]])

        texts = np.unique(np.array(texts), axis=0)
        text_dict = self.text_segregators(texts)

        for key, value in text_dict.items():
            t = ax.text(s=key, x=value[0], y=value[1], zorder=100)
            t.set_bbox(
                dict(
                    facecolor="white",
                    alpha=1.0,
                    pad=2.25,
                    edgecolor="black",
                    linewidth=0.5,
                )
            )

        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.set_ylim(-100, 4000)
        ax.set_xticks([])
        ax.set_ylabel("Temperature (K)", fontsize=12)
        ax.set_xlabel("Reaction Coordinate", fontsize=12)

        if self.save_flag:
            self.save_figure(
                folders=["pathways_plots"],
                file_name="-".join(sorted(self.composition)),
                fig=fig,
            )

        return ax, fig
