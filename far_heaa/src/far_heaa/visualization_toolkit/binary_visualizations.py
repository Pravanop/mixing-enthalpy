import numpy as np
from matplotlib import pyplot as plt

from typing import Tuple, Literal
from far_heaa.visualization_toolkit.visualizations import Visualizations


class binaryVizualization(Visualizations):
    """
    A class to visualize the binary phase diagram of a given alloy system. It calculates and plots
    the miscibility temperatures across a mole fraction grid for a binary composition.

    Args:
            composition (list[str]): A list containing the two elements in the binary alloy.
            lattice (str): The lattice type to consider (e.g., 'FCC', 'BCC').
            meta_data (dict): A dictionary containing metadata such as grid size, flags, and file paths.
            save_flag (bool): A flag indicating whether the plot should be saved.

    Raises:
            ValueError: If other than two elements are provided for the composition.
    Methods:
            find_misc_temperatures: Finds the miscibility temperatures for the binary alloy across a mole fraction grid.
            plot_misc_temperatures: Plots the miscibility temperatures for the binary alloy.
    """

    def __init__(
        self,
        composition: list[str],
        lattice: Literal["FCC", "BCC", "HCP", "min"],
        meta_data: dict,
        save_flag: bool,
    ):
        """
        Initialize the binaryVizualization class with the binary composition, lattice type, and metadata.

        Args:
                composition (list[str]): The two elements in the binary alloy.
                lattice (Literal['FCC', 'BCC', 'HCP', 'min']): The lattice structure. 'min' is whichever lattice gives minimum mixing enthalpy for that mole fraction.
                meta_data (dict): A dictionary containing metadata like grid size, file paths, and flags.
                save_flag (bool): A flag indicating whether the plot should be saved.

        Raises:
                ValueError: If the composition does not contain exactly two elements.
        """

        super().__init__(lattice, meta_data)

        self.composition = composition

        if len(self.composition) != 2:
            raise ValueError("Only provide binary compositions!")

        self.save_flag = save_flag

    def find_misc_temperatures(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Finds the miscibility temperatures for the binary alloy across a mole fraction grid.

        Returns:
                Tuple[np.ndarray, np.ndarray]: A tuple containing:
                        - mol_grid (np.ndarray): The mole fraction grid for the binary alloy.
                        - misc_temp (np.ndarray): The miscibility temperatures corresponding to each composition.

        Example::
        
                mol_grid, misc_temp = self.find_misc_temperatures()
        """
        mol_grid, misc_temp = self.grid_iterator.misc_temperature_across_grid(
            composition=self.composition,
            mol_grid_size=self.mol_grid_size,
            lattice=self.lattice,
            phase_flag=True,
        )
        return mol_grid[:-1], misc_temp[:-1]

    def plot_misc_temperatures(self) -> Tuple[plt.Axes, plt.Figure]:
        """
        Plots the miscibility temperatures for the binary alloy.

        Returns:
                Tuple[plt.Axes, plt.Figure]: A tuple containing:
                        - ax (plt.Axes): The matplotlib axes object for the plot.
                        - fig (plt.Figure): The matplotlib figure object for the plot.

        Example::
        
                ax, fig = self.plot_misc_temperatures()
        """
        mol_grid, misc_temp = self.find_misc_temperatures()

        # Replace -1 values with 5000 for visualization
        for idx, temp in enumerate(misc_temp):
            if temp == -1:
                misc_temp[idx] = 5000

        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        ax.set_facecolor("#004488")
        ax.plot(
            mol_grid[:, 0], misc_temp, color="black", linewidth=2, label="_nolegend_"
        )
        ax.fill(mol_grid[:, 0], misc_temp, color="#DDAA33", label="Stable")

        T_m = [0, 0]
        T_m[1] = self.tm.avg_T_melt(self.composition[0], [1])
        T_m[0] = self.tm.avg_T_melt(self.composition[1], [1])
        ax.plot([0, 1], T_m, "--", color="red", linewidth=1.5)

        ax.set_xlim(0, 1)
        ax.set_xticks([0, 1], sorted(self.composition, reverse=True))
        ax.set_ylabel("T (K)")
        ax.set_xlabel("x")
        ax.set_yticks(np.linspace(0, 3000, 6))
        ax.set_ylim(0, max(T_m) + 100)

        self.composition.reverse()
        ax.set_xticks([0, 1], self.composition)
        ax.legend(
            ["Unstable", "Average T$_{melt}$"],
            ncols=2,
            bbox_to_anchor=(1.1, 1.15),
            frameon=False,
        )

        # Adjust plot layout
        plt.subplots_adjust(left=0.22, right=0.95, bottom=0.15)

        # Save the plot if save_flag is True
        if self.save_flag:
            self.save_figure(
                folders=["binary_phase_diagram"],
                file_name="-".join(sorted(self.composition, reverse=True)),
                fig=fig,
            )

        return ax, fig
