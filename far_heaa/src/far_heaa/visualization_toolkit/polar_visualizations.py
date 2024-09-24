import numpy as np
from matplotlib import pyplot as plt, cm
from matplotlib.colors import Normalize
from typing import Tuple, List, Literal

from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.grids_and_combinations.grid_creation import CompositionGrid
from far_heaa.math_operations.polar_calculations import PolarMaths as pm
from far_heaa.data_processing.fancy_list_extractions import FancyListExtractions as fLE
from far_heaa.visualization_toolkit.visualizations import Visualizations


class PolarVisualizations(Visualizations):
    """
    A class for generating polar visualizations of miscibility temperatures and related properties for a given alloy system.

    Methods:
        misc_temp(member_pos, x, N, flag):
            Retrieves miscibility temperature data for the specified composition and lattice.

        draw_circle_in_polar(radius, ax):
            Draws a circle with the given radius on a polar plot.

        get_n_colors_from_cmap(cmap_name, N):
            Retrieves N equally spaced colors from the specified colormap.

        text_flipper(angle):
            Flips text orientation based on the angle in a polar plot.

        scatter_center(scatter, ax):
            Plots the central scatter point in a polar plot.

        plot_line(angle_degrees, x_values, temp_list, ax, cmap, norm, zorder):
            Plots a line with varying color on a polar plot based on temperature values.

        plot_colored_secant(ax, radius, angle_start, angle_end, misc_temp_list_sec, cmap, norm, zorder=50):
            Plots a colored secant line between two angles on a polar plot.

        set_ax_params(ax):
            Sets up the polar plot with proper parameters like colorbars and labels.

        make_one_bar(ax, n_alloy, member_pos, i, N, angle, idx2):
            Plots a single bar representing miscibility temperature on a polar plot.

        total_plot():
            Generates a full polar plot for all compositions in the alloy system.

        plot_N(N_ind, transmute_indices):
            Plots polar diagrams for a specific N-index with optional transmutation lines.
    """

    def __init__(
        self,
        composition: List[str],
        lattice: Literal["FCC", "BCC", "HCP", "min"],
        meta_data: dict,
        save_flag: bool,
        type_flag: Literal["misc_T", "e_hull"],
        **kwargs,
    ):
        """
        Initializes the PolarVisualizations class with composition, lattice, metadata, and save flag.

        Args:
            composition (List[str]): The alloy composition.
            lattice (Literal['FCC', 'BCC', 'HCP', 'min']): The lattice structure. 'min' is whichever lattice gives minimum mixing enthalpy for that mole fraction.
            meta_data (dict): Metadata containing grid size, folder paths, etc.
            save_flag (bool): A flag indicating whether to save the visualizations.
        """
        super().__init__(lattice=lattice, meta_data=meta_data)
        self.composition = composition
        if len(self.composition) < 3:
            raise ValueError(
                "Composition should have at least 3 elements. Maybe you are better off using other visualizations."
            )

        self.cmap = cm.get_cmap("plasma")
        self.save_flag = save_flag

        self.type_flag = type_flag
        if self.type_flag == "misc_T":
            t_max = max([self.tm.avg_T_melt(i, mol_ratio=[]) for i in self.composition])
            self.temp_gradation = kwargs.get("temp_gradation", 25)
            self.temp_grid = list(np.linspace(0, t_max, self.temp_gradation))
            self.conv_hull = self.grid_iterator.temp_iterator(
                composition=self.composition, temp_grid=self.temp_grid
            )
            self.norm = Normalize(vmin=0, vmax=t_max)

        elif self.type_flag == "e_hull":
            self.norm = Normalize(vmin=0, vmax=0.25)

        self.mol_gradation = kwargs.get("mol_gradation", 15)
        self.x = np.linspace(0, 1, self.mol_gradation)
        self.figsize = (len(self.composition) + 1, len(self.composition) + 1)
        self.linewidth = 9.5 - (len(self.composition) - 3)
        self.y_bias = 0.4

    def misc_temp(
        self,
        member_pos: List[int],
        x: np.ndarray,
        N: int,
        flag: Literal["add", "transmutate"],
    ) -> np.ndarray:
        """
        Retrieves miscibility temperature data for the specified composition and lattice.

        Args:
            member_pos (List[int]): Indices of the members in the composition.
            x (np.ndarray): The mole fraction grid.
            N (int): The dimensionality of the composition.
            flag (Literal['add', 'transmutate']): A flag to determine whether to "add" or "transmute" elements.

        Returns:
            np.ndarray: Miscibility temperatures for the given composition.
        """
        n = len(self.composition)

        mol_grid = (
            CompositionGrid.create_high_sym_mol_grid(
                x=x, n=n, N=N, change_idx=member_pos
            )
            if flag == "add"
            else CompositionGrid.create_mol_grid_transmutation(
                x=x, n=n, transmutation_indice=member_pos
            )
        )

        mol_grid, misc_temp = self.grid_iterator.misc_temperature_across_grid(
            composition=self.composition,
            mol_grid_size=mol_grid,
            lattice=self.lattice,
            phase_flag=False,
            conv_hull=self.conv_hull,
            temp_grid=self.temp_grid,
        )
        return misc_temp

    def e_hull(
        self,
        temperature: float,
        member_pos: List[int],
        x: np.ndarray,
        N: int,
        flag: Literal["add", "transmutate"],
    ):

        n = len(self.composition)

        mol_grid = (
            CompositionGrid.create_high_sym_mol_grid(
                x=x, n=n, N=N, change_idx=member_pos
            )
            if flag == "add"
            else CompositionGrid.create_mol_grid_transmutation(
                x=x, n=n, transmutation_indice=member_pos
            )
        )

        mol_grid, e_hull, _ = self.grid_iterator.e_hull_across_grid(
            composition=self.composition,
            mol_grid_size=mol_grid,
            lattice=self.lattice,
            single_temp_flag=True,
            temp_gradation=temperature,
        )
        return e_hull

    def draw_circle_in_polar(self, radius: float, ax: plt.Axes) -> None:
        """
        Draws a circle with the given radius on a polar plot.

        Args:
            radius (float): The radius of the circle.
            ax (plt.Axes): The polar plot axes on which the circle is drawn.
        """
        theta = np.linspace(0, 2 * np.pi, 100)
        ax.plot(
            theta,
            [radius + self.y_bias] * len(theta),
            linewidth=0.8,
            zorder=0,
            color="black",
            linestyle="--",
            alpha=0.1,
        )

    @staticmethod
    def get_n_colors_from_cmap(cmap_name: str, N: int) -> List:
        """
        Retrieves N equally spaced colors from the specified colormap.

        Args:
            cmap_name (str): The name of the colormap.
            N (int): The number of colors to retrieve.

        Returns:
            List: A list of colors from the colormap.
        """
        cmap = plt.get_cmap(cmap_name)
        return [cmap(i) for i in np.linspace(0, 1, N)]

    @staticmethod
    def text_flipper(angle: float) -> float:
        """
        Flips text orientation based on the angle in a polar plot.

        Args:
            angle (float): The angle in degrees.

        Returns:
            float: The adjusted angle for text orientation.
        """
        return angle + 180 if 90 < angle < 270 else angle

    def scatter_center(self, scatter: float, ax: plt.Axes) -> None:
        """
        Plots the central scatter point in a polar plot.

        Args:
            scatter (float): The value to normalize the color of the scatter point.
            ax (plt.Axes): The polar plot axes.
        """
        
        theta = np.linspace(0, 2 * np.pi, 100)
        # self.draw_circle_in_polar(radius=self.y_bias, ax=ax)
        ax.plot(
            theta,
            [self.y_bias - 0.05] * len(theta),
            linewidth=1.5,
            zorder=100,
            color="black",
            linestyle="-",
            alpha=1,
        )
        ax.fill(
            theta,
            [self.y_bias - 0.05] * len(theta),
            zorder=100,
            color=self.cmap(self.norm(scatter)),
            alpha=1,
        )

    def plot_line(
        self,
        angle_degrees: float,
        x_values: np.ndarray,
        temp_list: np.ndarray,
        ax: plt.Axes,
        cmap,
        norm,
        zorder: int,
    ) -> None:
        """
        Plots a line with varying color on a polar plot based on temperature values.

        Args:
            angle_degrees (float): The angle in degrees where the line should be drawn.
            x_values (np.ndarray): The x-values corresponding to mole fractions.
            temp_list (np.ndarray): List of temperature values for coloring the line.
            ax (plt.Axes): The polar plot axes.
            cmap: Colormap to use for coloring the line.
            norm: Normalization for the colormap.
            zorder (int): The drawing order.
        """

        angle_radians = np.radians(angle_degrees)
        for i in range(len(x_values)):
            ax.vlines(
                x=angle_radians,
                ymax=x_values[i] + self.y_bias,
                ymin=self.y_bias,
                color=cmap(norm(temp_list[i])),
                linewidth=self.linewidth,
                zorder=zorder + len(x_values) - i,
            )

    def plot_colored_secant(
        self,
        ax: plt.Axes,
        radius: float,
        angle_start: float,
        angle_end: float,
        misc_temp_list_sec: np.ndarray,
        cmap,
        norm,
        zorder: int = 50,
    ) -> None:
        """
        Plots a colored secant line between two angles on a polar plot.

        Args:
            ax (plt.Axes): The polar plot axes.
            radius (float): The radius for the secant line.
            angle_start (float): The starting angle in radians.
            angle_end (float): The ending angle in radians.
            misc_temp_list_sec (np.ndarray): The miscibility temperature values for each segment of the secant.
            cmap: Colormap for coloring the secant.
            norm: Normalization for the colormap.
            zorder (int, optional): The drawing order. Defaults to 50.
        """
        n_segments = len(misc_temp_list_sec)
        x_start, y_start = pm.polar_to_cartesian(radius, angle_start)
        x_end, y_end = pm.polar_to_cartesian(radius, angle_end)
        total_x = np.linspace(x_start, x_end, n_segments)
        slope = (y_end - y_start) / (x_end - x_start)
        y = slope * (total_x - x_start) + y_start
        radius, theta = pm.cartesian_to_polar(total_x, y)

        for idx in range(len(misc_temp_list_sec)):
            if idx + 1 < len(misc_temp_list_sec):
                x_coord = [theta[idx], theta[idx + 1]]
                y_coord = [radius[idx], radius[idx + 1]]
                ax.plot(
                    x_coord,
                    y_coord,
                    c=cmap(norm(misc_temp_list_sec[idx])),
                    linewidth=self.linewidth,
                    zorder=zorder + len(misc_temp_list_sec) - idx,
                )

    def set_ax_params(self, ax: plt.Axes) -> None:
        """
        Sets up the polar plot with proper parameters like colorbars and labels.

        Args:
            ax (plt.Axes): The polar plot axes.
        """
        ax.set_yticks([])
        ax.set_xticks([])
        ax.spines["polar"].set_visible(False)
        ax.grid(False)
        sm = plt.cm.ScalarMappable(cmap=self.cmap, norm=self.norm)
        sm.set_array([])  # We need this for colorbar to work
        cbar = plt.colorbar(
            sm, ax=ax, aspect=30, fraction=0.05, orientation="horizontal"
        )  # Unified colorbar
        if self.type_flag == "misc_T":
            cbar.set_label("$T_{misc}$ (K)", fontsize=12)
        else:
            cbar.set_label("$E_{hull}$ (eV/atom)", fontsize=12)

    def find_misc_temperatures(
        self, member_pos: List[int], N: int, flag: Literal["add", "transmutate"]
    ) -> np.ndarray:
        misc_temp_list = self.misc_temp(member_pos=member_pos, x=self.x, flag=flag, N=N)
        misc_temp_list = np.array([5000 if i == -1 else i for i in misc_temp_list])
        return misc_temp_list

    def make_one_bar(
        self,
        ax: plt.Axes,
        n_alloy: int,
        member_pos: List[int],
        i: str,
        N: int,
        angle: float,
        idx2: int,
        **kwargs,
    ) -> float:
        """
        Plots a single bar representing miscibility temperature on a polar plot.

        Args:
            ax (plt.Axes): The polar plot axes.
            n_alloy (int): The number of elements in the alloy.
            member_pos (List[int]): The position of the elements in the composition.
            i (str): The alloy composition string.
            N (int): The dimensionality of the alloy.
            angle (float): The angle in degrees for the bar placement.
            idx2 (int): Index for coloring the bar.

        Returns:
            float: The first miscibility temperature value for this bar.
        """
        line_colors = self.get_n_colors_from_cmap("Dark2", len(self.composition) - 1)
        temperature = kwargs.get("temperature", None)
        if self.type_flag == "e_hull" and not temperature:
            raise ValueError("Temperature should be provided for e_hull plot.")
        if temperature:
            misc_temp_list = self.e_hull(
                temperature=temperature,
                member_pos=member_pos,
                x=self.x,
                N=N,
                flag="add",
            )
        else:
            misc_temp_list = self.find_misc_temperatures(member_pos, N, flag="add")
        self.plot_line(
            angle_degrees=angle,
            x_values=self.x * pm.distance_calculator(n_alloy, N),
            temp_list=misc_temp_list,
            ax=ax,
            cmap=self.cmap,
            norm=self.norm,
            zorder=1,
        )
        angle_radians = np.radians(angle)
        ax.vlines(
            angle_radians,
            ymin=self.y_bias,
            ymax=(self.x * pm.distance_calculator(n_alloy, n_alloy - 1))[-1]
            + self.y_bias,
            linestyles="-",
            color=line_colors[idx2],
            zorder=0,
            alpha=0.7,
            linewidth=2,
        )
        rotation = self.text_flipper(angle=angle)
        ax.text(
            angle_radians,
            pm.distance_calculator(n_alloy, n_alloy - 1)
            + N * 0.2
            - 0.01 * n_alloy
            + self.y_bias,
            i,
            ha="center",
            va="center",
            color="black",
            rotation=rotation,
        )
        return float(misc_temp_list[0])

    def plot_total(self, **kwargs) -> Tuple[plt.Axes, plt.Figure]:
        """
        Generates a full polar plot for all compositions in the alloy system.

        Returns:
            Tuple[plt.Axes, plt.Figure]: The matplotlib axes and figure objects for the plot.
        """
        n_alloy = len(self.composition)
        fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=self.figsize)

        angles = pm.divide_circle_degrees(pm.total_num_bars(n_alloy))
        count = 0
        scatter = 0
        temperature = kwargs.get("temperature", None)
        print(temperature)
        if self.type_flag == "e_hull" and temperature is None:
            raise ValueError("Temperature should be provided for e_hull plot.")
        for idx2, N in enumerate(range(1, len(self.composition))):

            combs = (
                self.composition
                if N == 1
                else list(
                    MultinaryCombinations.create_multinary(
                        element_list=self.composition, no_comb=[N]
                    ).values()
                )[0]
            )

            for idx, i in enumerate(combs):
                angle = angles[count]
                temp_i = i.split("-")
                member_pos = fLE.find_indices(self.composition, temp_i)
                if temperature:
                    scatter = self.make_one_bar(
                        ax,
                        n_alloy,
                        member_pos,
                        i,
                        N,
                        float(angle),
                        idx2,
                        temperature=temperature,
                    )
                else:
                    scatter = self.make_one_bar(
                        ax, n_alloy, member_pos, i, N, float(angle), idx2
                    )

                count += 1

        for N in range(1, len(self.composition)):
            self.draw_circle_in_polar(radius=pm.distance_calculator(n_alloy, N), ax=ax)

        self.scatter_center(scatter, ax)
        self.set_ax_params(ax)

        if self.save_flag:
            self.save_figure(
                folders=["polar_plots", self.type_flag, "total"],
                file_name=f'{"".join(sorted(self.composition))}',
                fig=fig,
            )

        return ax, fig

    def plot_subset(
        self, N_ind: int, transmute_indices: List[int], **kwargs
    ) -> Tuple[plt.Axes, plt.Figure]:
        """
        Plots polar diagrams for a specific N-index with optional transmutation lines.

        Args:
            N_ind (int): The N-index to plot.
            transmute_indices (List[int]): The indices of elements to transmute.

        Returns:
            Tuple[plt.Axes, plt.Figure]: The matplotlib axes and figure objects for the plot.
        """
        n_alloy = len(self.composition)
        fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=self.figsize)

        melt_T = []
        for idx, i in enumerate(self.composition):
            mol_ratio = [1 if idx == i else 0 for i in range(n_alloy)]
            melt_T.append(self.tm.avg_T_melt(composition=[i], mol_ratio=mol_ratio))

        count = 0
        if N_ind == 1:
            combs = self.composition
        else:
            combs = list(
                MultinaryCombinations.create_multinary(
                    element_list=self.composition, no_comb=[N_ind], sort=False
                ).values()
            )[0]

        combs_n_N = [
            "-".join(list(set(self.composition).difference(set(i.split("-")))))
            for i in combs
        ]

        scatter = 0
        temperature = kwargs.get("temperature", None)
        if self.type_flag == "e_hull" and not temperature:
            raise ValueError("Temperature should be provided for e_hull plot.")
        for idx2, comb in enumerate([combs, combs_n_N]):
            angles = pm().angle_assigner(length=len(comb))

            for idx, i in enumerate(comb):
                angle = angles[idx] + idx2 * 180
                temp_i = i.split("-")
                member_pos = fLE.find_indices(self.composition, temp_i)
                N = len(temp_i)
                if temperature:
                    scatter = self.make_one_bar(
                        ax,
                        n_alloy,
                        member_pos,
                        i,
                        N,
                        float(angle),
                        idx2,
                        temperature=temperature,
                    )
                else:
                    scatter = self.make_one_bar(
                        ax, n_alloy, member_pos, i, N, float(angle), idx2
                    )

                angle_radians = np.radians(angle)
                if transmute_indices and self.type_flag == "misc_T":
                    if idx2 != 0 and idx == min(transmute_indices):
                        misc_temp_list_sec = self.find_misc_temperatures(
                            transmute_indices, N, flag="transmutate"
                        )
                        self.plot_colored_secant(
                            ax,
                            float(
                                (self.x * pm.distance_calculator(n_alloy, n_alloy - 1))[
                                    -1
                                ]
                                + self.y_bias
                            ),
                            angle_radians,
                            np.radians(angles[transmute_indices[1]] + idx2 * 180),
                            misc_temp_list_sec,
                            cmap=self.cmap,
                            norm=self.norm,
                        )
                self.draw_circle_in_polar(
                    radius=pm.distance_calculator(n_alloy, N), ax=ax
                )
                count += 1

        self.scatter_center(scatter, ax)
        self.set_ax_params(ax)

        if self.save_flag:
            self.save_figure(
                folders=["polar_plots", self.type_flag, f"{N_ind}"],
                file_name=f'{"".join(sorted(self.composition))}',
                fig=fig,
            )

        return ax, fig
