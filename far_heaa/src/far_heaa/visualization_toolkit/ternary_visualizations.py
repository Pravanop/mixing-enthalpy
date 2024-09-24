import numpy as np
from matplotlib import pyplot as plt
import mpltern
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from typing import Tuple, Literal

from far_heaa.math_operations.thermo_calculations import ThermoMaths
from far_heaa.phase_diagram.grid_iterators import GridIterator
from far_heaa.visualization_toolkit.visualizations import Visualizations


class TernaryVisualization(Visualizations):
    """
    A class for visualizing ternary phase diagrams and miscibility temperatures for a given alloy composition.

    Args:
            composition (list[str]): A list of three elements in the ternary alloy system.
            lattice (Literal['FCC', 'BCC', 'HCP', 'min']): The lattice structure. 'min' is whichever lattice gives minimum mixing enthalpy for that mole fraction.
            meta_data (dict): A dictionary containing metadata such as grid size, file paths, and flags.
            save_flag (bool): A flag indicating whether to save the visualizations.
            contour_flag (bool): A flag indicating whether to use contour plots or scatter plots.

    Attributes:
            cmap (str): The colormap to use for plotting.
            tm (ThermoMaths): An instance of the ThermoMaths class for thermodynamic calculations.
            grid_iterator (GridIterator): An instance of the GridIterator class for grid-based iteration.
            save_flag (bool): A flag indicating whether to save the generated plots.
            contour_flag (bool): A flag indicating whether to generate contour plots.
    Methods:
            find_misc_temperatures: Finds the miscibility temperatures for the ternary alloy.
            find_isotherm: Finds the isotherm (temperature slice) of the energy above the convex hull for the ternary alloy.
            plot_isotherm: Plots the isotherm (temperature slice) of the energy above the convex hull for the ternary alloy.
            plot_misc_temperatures: Plots the miscibility temperatures for the ternary alloy.
            plot_ternary_visualizations: Generates and saves ternary visualizations, including miscibility temperature plots and isotherms for a given temperature range and gradation.
    """

    def __init__(
        self,
        composition: list[str],
        lattice: Literal["FCC", "BCC", "HCP", "min"],
        meta_data: dict,
        save_flag: bool,
        contour_flag: bool,
    ):
        """
        Initializes the TernaryVizualization class with the ternary alloy composition, lattice type, and metadata.

        Args:
                composition (list[str]): A list of three elements in the ternary alloy system.
                lattice (Literal['FCC', 'BCC', 'HCP', 'min']): The lattice structure. 'min' is whichever lattice gives minimum mixing enthalpy for that mole fraction.
                meta_data (dict): A dictionary containing metadata such as grid size, file paths, and flags.
                save_flag (bool): A flag indicating whether to save the visualizations.
                contour_flag (bool): A flag indicating whether to use contour plots or scatter plots.

        Raises:
                ValueError: If the provided composition does not contain exactly three elements.
        """
        super().__init__(lattice=lattice, meta_data=meta_data)

        self.cmap = "plasma"

        self.composition = composition
        if len(self.composition) != 3:
            raise ValueError("Only provide ternary compositions!")

        self.save_flag = save_flag
        self.contour_flag = contour_flag

    def find_misc_temperatures(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Finds miscibility temperatures for the ternary alloy across a mole fraction grid.

        Returns:
                Tuple[np.ndarray, np.ndarray]: A tuple containing:
                        - mol_grid (np.ndarray): The mole fraction grid for the ternary alloy.
                        - misc_temp (np.ndarray): The miscibility temperatures corresponding to each composition.
        """
        mol_grid, misc_temp = self.grid_iterator.misc_temperature_across_grid(
            composition=self.composition,
            mol_grid_size=self.mol_grid_size,
            lattice=self.lattice,
            phase_flag=True,
        )

        return mol_grid, misc_temp

    def find_isotherm(self, temperature: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Finds the isotherm (temperature slice) of the energy above the convex hull for the ternary alloy.

        Args:
                temperature (float): The temperature in Kelvin to evaluate the isotherm.

        Returns:
                Tuple[np.ndarray, np.ndarray]: A tuple containing:
                        - mol_grid (np.ndarray): The mole fraction grid for the ternary alloy.
                        - stables (np.ndarray): The stability of the alloy compositions at the given temperature.
        """
        mol_grid, stables, _ = self.grid_iterator.e_hull_across_grid(
            composition=self.composition,
            mol_grid_size=self.mol_grid_size,
            lattice=self.lattice,
            single_temp_flag=True,
            temp_gradation=temperature,
        )
        return mol_grid, stables

    def plot_isotherm(self, temperature: float) -> Tuple[plt.Axes, plt.Figure]:
        """
        Plots the isotherm (temperature slice) of the energy above the convex hull for the ternary alloy.

        Args:
                temperature (float): The temperature in Kelvin to evaluate the isotherm.

        Returns:
                Tuple[plt.Axes, plt.Figure]: A tuple containing the matplotlib axes and figure objects for the plot.
        """
        mol_grid, stables = self.find_isotherm(temperature)

        for idx, stable in enumerate(stables):
            if np.isclose(stable, 0.0, atol=1e-3):
                stables[idx] = 0
            else:
                stables[idx] = 1

        t, l, r = mol_grid[:, 0], mol_grid[:, 1], mol_grid[:, 2]
        fig = plt.figure()

        ax = fig.add_subplot(projection="ternary")
        ax.grid()
        norm = Normalize(vmin=0, vmax=0.25, clip=False)
        if not self.contour_flag:
            ax.scatter(t, l, r, c=stables, cmap=self.cmap, marker="h", s=60, norm=norm)
        else:
            data = np.concatenate(
                [
                    t.reshape(-1, 1),
                    l.reshape(-1, 1),
                    r.reshape(-1, 1),
                    stables.reshape(-1, 1),
                ],
                axis=1,
            )
            ax.tricontourf(
                data[:, 0],
                data[:, 1],
                data[:, 2],
                data[:, 3],
                cmap=self.cmap,
                norm=norm,
            )

        ax.grid(False)
        ax.set_tlabel(f"{self.composition[0]}")
        ax.set_llabel(f"{self.composition[1]}")
        ax.set_rlabel(f"{self.composition[2]}")

        if self.save_flag:
            if self.contour_flag:
                self.save_figure(
                    folders=[
                        "ternary_phase_diagram",
                        "isotherms",
                        "contours",
                        f"{'-'.join(sorted(self.composition, reverse=True))}",
                    ],
                    file_name=f"{temperature}",
                    fig=fig,
                )
            else:
                self.save_figure(
                    folders=[
                        "ternary_phase_diagram",
                        "isotherms",
                        "scatters",
                        f"{'-'.join(sorted(self.composition, reverse=True))}",
                    ],
                    file_name=f"{temperature}",
                    fig=fig,
                )

        return ax, fig

    def plot_misc_temperatures(self) -> Tuple[plt.Axes, plt.Figure]:
        """
        Plots the miscibility temperatures for the ternary alloy.

        Returns:
                Tuple[plt.Axes, plt.Figure]: A tuple containing the matplotlib axes and figure objects for the plot.
        """
        mol_grid, misc_temp = self.find_misc_temperatures()

        for idx, temp in enumerate(misc_temp):
            if temp == -1:
                misc_temp[idx] = 5000

        t, l, r = mol_grid[:, 0], mol_grid[:, 1], mol_grid[:, 2]
        fig = plt.figure()

        ax = fig.add_subplot(projection="ternary")
        ax.grid()
        cax = ax.inset_axes((1.03, 0.1, 0.05, 0.9), transform=ax.transAxes)
        norm = Normalize(vmin=0, vmax=3600, clip=False)
        if not self.contour_flag:
            ax.scatter(t, l, r, c=misc_temp, cmap=self.cmap, marker="h", s=60)
        else:
            data = np.concatenate(
                [
                    t.reshape(-1, 1),
                    l.reshape(-1, 1),
                    r.reshape(-1, 1),
                    misc_temp.reshape(-1, 1),
                ],
                axis=1,
            )
            ax.tricontourf(
                data[:, 0],
                data[:, 1],
                data[:, 2],
                data[:, 3],
                cmap=self.cmap,
                norm=norm,
            )

        sm = ScalarMappable(cmap=self.cmap, norm=norm)
        colorbar = fig.colorbar(sm, cax=cax)
        colorbar.set_label("$T_{misc}$ K", rotation=270, va="baseline")
        ax.grid(False)
        ax.set_tlabel(f"{self.composition[0]}")
        ax.set_llabel(f"{self.composition[1]}")
        ax.set_rlabel(f"{self.composition[2]}")

        if self.save_flag:
            if self.contour_flag:
                self.save_figure(
                    folders=["ternary_phase_diagram", "misc_temp", "contours"],
                    file_name=f"{'-'.join(sorted(self.composition, reverse=True))}",
                    fig=fig,
                )

            else:
                self.save_figure(
                    folders=["ternary_phase_diagram", "misc_temp", "scatter"],
                    file_name=f"{'-'.join(sorted(self.composition, reverse=True))}",
                    fig=fig,
                )

        return ax, fig

    def plot_ternary_visualizations(
        self,
        T_min: float,
        T_max: float,
        T_gradation: float,
        switch_contour_off: bool = False,
        switch_scatter_off: bool = False,
    ) -> None:
        """
        Generates and saves ternary visualizations, including miscibility temperature plots and isotherms
        for a given temperature range and gradation. (both contours and scatters

        Args:
                T_min (float): The minimum temperature for the isotherm grid.
                T_max (float): The maximum temperature for the isotherm grid.
                T_gradation (float): The temperature step/gradation between isotherms.
                switch_contour_off (bool): Whether to switch the contour off. Default is False.
                switch_scatter_off (bool): Whether to switch the scatter off. Default is False.

        RaisesError: If switch_contour_off and switch_scatter_off is True.

        Example::
        
                self.plot_ternary_visualizations(T_min=300, T_max=1500, T_gradation=100)
        """
        # Create a grid of temperatures for isotherms
        isotherm_grid = np.arange(T_min, T_max + 100, T_gradation)

        # Set flags for saving and for scatter plots
        self.save_flag = True

        if switch_scatter_off and switch_contour_off:
            raise "Why do you want to waste your time with this function?"

        if not switch_scatter_off:
            self.contour_flag = False

            # Plot miscibility temperatures with scatter plots
            _, _ = self.plot_misc_temperatures()

            # Plot each isotherm using scatter plot
            for isotherm in isotherm_grid:
                _, _ = self.plot_isotherm(temperature=isotherm)

        if not switch_contour_off:
            # Switch to contour plots
            self.contour_flag = True

            # Plot miscibility temperatures with contour plots
            _, _ = self.plot_misc_temperatures()

            # Plot each isotherm using contour plot
            for isotherm in isotherm_grid:
                _, _ = self.plot_isotherm(temperature=isotherm)
