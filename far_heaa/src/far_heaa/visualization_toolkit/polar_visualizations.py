import pickle

import numpy as np
from matplotlib import pyplot as plt, cm
from matplotlib.colors import Normalize
from typing import Tuple, List, Literal

from matplotlib.patches import Rectangle, Wedge
import matplotlib.colors as mcolors
from matplotlib.pyplot import get_cmap
from pybtex.textutils import width
from pymatgen.core import Element

from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.grids_and_combinations.grid_creation import CompositionGrid
from far_heaa.math_operations.polar_calculations import PolarMaths as pm
from far_heaa.data_processing.fancy_list_extractions import FancyListExtractions as fLE
from far_heaa.visualization_toolkit.visualizations import Visualizations
import matplotlib

# matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.family'] = 'Helvetica'
matplotlib.rcParams['font.size'] = 10

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
        type_flag: Literal["misc_T", "e_hull", "entropy", "melt"],
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
        self.save_flag = save_flag
        self.fontsize = 10
        self.type_flag = type_flag
        self.mol_dict = []
        if self.type_flag == "misc_T":
            n = 10
            base_cmap = plt.get_cmap('plasma_r')
            colors = [
                "darkgrey",
                *base_cmap(np.linspace(0, 1, n))
            ]
            custom_cmap = mcolors.ListedColormap(colors)
            self.cmap = custom_cmap
            boundaries = [-3000, 200] + list(np.round(np.linspace(200, 3750, n)))
            self.norm = mcolors.BoundaryNorm(boundaries, self.cmap.N)


            t_max = self.tm.avg_T_melt(self.composition, mol_ratio=[1/len(self.composition)]*len(self.composition))
            self.temp_gradation = kwargs.get("temp_gradation", 25)
            self.temp_grid = list(np.linspace(0, t_max, self.temp_gradation))
            self.conv_hull = self.grid_iterator.temp_iterator(
                composition=self.composition, temp_grid=self.temp_grid
            )


        elif self.type_flag == "e_hull":
            base_cmap = plt.get_cmap('OrRd')
            n = 10
            colors = [
                'cornflowerblue',  # Dark blue for 0
                *base_cmap(np.linspace(0.2, 0.6, n)),  # OrRd for values between 0.0 (exclusive) and 0.05
                'darkred'  # Darkest red for values above 0.05
            ]

            # Create a custom colormap
            custom_cmap = mcolors.ListedColormap(colors)

            # Define boundaries
            boundaries = [0.0, 0.005] + list(np.round(np.linspace(0.0051, 0.05, n), 3)) + [0.4]

            # Create a BoundaryNorm for normalization
            norm = mcolors.BoundaryNorm(boundaries, custom_cmap.N)

            # Assign colormap and norm to your object
            self.cmap = custom_cmap
            self.norm = norm

        elif self.type_flag == "gibbs":
            base_cmap = plt.get_cmap('coolwarm')
            n = 10
            colors = [
                *base_cmap(np.linspace(0, 1, n))
            ]
            custom_cmap = mcolors.ListedColormap(colors)
            self.cmap = custom_cmap
            boundaries = list(np.round(np.linspace(-50, 50, n),0))
            self.norm = mcolors.BoundaryNorm(boundaries, self.cmap.N)

        elif self.type_flag == "entropy":
            self.cmap = plt.get_cmap('cubehelix', 8)
            self.cmap = mcolors.ListedColormap(self.cmap(np.linspace(0.2, 1, 7))[:-1])
            limit = -(1/len(self.composition))*np.log(1/len(self.composition))*len(self.composition)
            self.norm = Normalize(vmin=0, vmax=limit)

        elif self.type_flag == "melt":
            t_min = 5000
            for j in self.composition:
                t_min_temp = self.tm.avg_T_melt(j, mol_ratio=[1/len(self.composition)]*len(self.composition))
                if t_min_temp < t_min:
                    t_min = t_min_temp
            self.cmap = plt.get_cmap('plasma')
            self.norm = Normalize(vmin=t_min, vmax=3700)
        elif self.type_flag == "density":
            self.cmap = plt.get_cmap('jet', 8)
            self.cmap = mcolors.ListedColormap(self.cmap(np.linspace(0.2, 1, 7))[:-1])
            self.norm = Normalize(vmin = 7, vmax= 9)
        elif self.type_flag == "elastic":
            self.cmap = plt.get_cmap('jet')
            self.norm = Normalize(vmin = 100, vmax = 300)

        self.mol_gradation = kwargs.get("mol_gradation", 15)
        self.x = np.linspace(0, 1, self.mol_gradation)
        # self.figsize = (len(self.composition) + 1, len(self.composition) + 1)
        self.figsize = (6.34, 3.33)
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
            is_differential=True
        )
        return misc_temp


    def elastic(self,
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

        elastic = []
        for j in mol_grid:
            modulus = np.array([Element(i).youngs_modulus for i in self.composition])
            elastic.append(np.dot(modulus, np.array(j)))
        return elastic

    def density(self,
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

        density = []
        for j in mol_grid:
            modulus = np.array([Element(i).atomic_mass/Element(i).molar_volume for i in self.composition])
            density.append(np.dot(modulus, np.array(j)))
        return density

    def melt_T(self,
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

        melt_T = []
        for i in mol_grid:
            melt_T.append(
                self.grid_iterator.tm.avg_T_melt(composition = self.composition,
                                   mol_ratio= i)
            )
        return melt_T

    def entropy(self,
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

        entropy = []
        for i in mol_grid:
            entropy.append(
                self.grid_iterator.tm.calc_config_entropy(mol_ratio=dict(zip(self.composition, list(i + 1e-8))))/self.tm.kb
            )

        print(entropy)
        return entropy

    def gibbs(self,
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

        gibbs = []
        omega1, omega2, omega3, omega4 = 0.420, 0.420, 0.420, 0.420
        for i in mol_grid:
            entropy = self.grid_iterator.tm.calc_config_entropy(mol_ratio=dict(zip(self.composition, list(i+1e-9))))
            hmix = sum([omega1*i[0]*i[1], omega2/2*i[0]*i[2], omega3/4*i[0]*i[3], omega4*i[1]*i[2], omega1*i[1]*i[3], omega3*i[2]*i[3]])
            gibbs.append(np.round((hmix - 1000*entropy)*1000, 0))

        print(mol_grid[-1], gibbs)
        return gibbs

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
        self.mol_dict.append(mol_grid)
        # mol_grid, e_hull, _ = self.grid_iterator.e_hull_across_grid(
        #     composition=self.composition,
        #     mol_grid_size=mol_grid,
        #     lattice=self.lattice,
        #     single_temp_flag=True,
        #     temp_gradation=temperature,
        # )
        return mol_grid

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

    #deprecated
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

    def plot_rectangle(
            self,
            angle_degrees: float,
            x_values: np.ndarray,
            temp_list: np.ndarray,
            ax: plt.Axes,
            cmap,
            norm,
            zorder: int,
            bar_width: float
    ) -> None:
        """
        Plots rectangles with varying color on a polar plot based on temperature values.

        Args:
            angle_degrees (float): The angle in degrees where the rectangles should be drawn.
            x_values (np.ndarray): The x-values corresponding to mole fractions.
            temp_list (np.ndarray): List of temperature values for coloring the rectangles.
            ax (plt.Axes): The polar plot axes.
            cmap: Colormap to use for coloring the rectangles.
            norm: Normalization for the colormap.
            zorder (int): The drawing order.
            bar_width (float): The width of the bars.
        """

        # angle_radians = np.round(np.radians(angle_degrees), decimals=2)
        # bar_width = np.radians(10) - np.radians(len(self.composition)*1.1)  # Adjust width as needed, here it's 1 degree
        #
        # for i in range(len(x_values)):
        #     # Define the position and dimensions of the rectangle
        #     rect = Rectangle(
        #         (angle_radians - bar_width / 2, self.y_bias),  # (theta, r) starting point
        #         width=bar_width,
        #         height=x_values[i],  # Height is based on the x_values
        #         facecolor=cmap(norm(temp_list[i])),
        #         zorder=zorder + len(x_values) - i,
        #         edgecolor=(0, 0, 0, 0.5),
        #         linewidth=0.5,
        #         antialiased=True
        #     )
        #     ax.add_patch(rect)
        angle_radians = np.radians(angle_degrees)
        
        bar_width = np.radians(360/bar_width - 2)
        for i in range(len(x_values) - 1):
            # Define the start and end of each radial segment
            if i == 0:
                inner_radius = x_values[i] + self.y_bias
            else:
                inner_radius = x_values[i] + self.y_bias
                
            outer_radius = x_values[i + 1] + self.y_bias
            
            # Create a rectangle to represent each segment in the bar
            rect = Rectangle(
                (angle_radians - bar_width/2, inner_radius),  # (theta, r) starting point
                width=bar_width,  # Width in radians (angle)
                height=outer_radius - inner_radius,  # Radial height of each segment
                facecolor=cmap(norm(temp_list[i+1])),
                edgecolor="#555555",
                linewidth=0.5,# No edge for a smooth gradient effect
                zorder=zorder
            )
            
            # Add the rectangle to the polar plot
            ax.add_patch(rect)

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

    def plot_colored_secant_rectangle(
            self,
            ax: plt.Axes,
            radius: float,
            angle_start: float,
            angle_end: float,
            misc_temp_list_sec: np.ndarray,
            cmap,
            norm,
            rect_height: float,
            zorder: int = 20,
    ) -> None:
        """
         Plots colored rectangles along an arc of a circle's circumference.

    Args:
        ax (plt.Axes): The polar plot axes.
        radius (float): The radius for placing the rectangles along the arc.
        angle_start (float): The starting angle in radians.
        angle_end (float): The ending angle in radians.
        misc_temp_list_sec (np.ndarray): The miscibility temperature values for each rectangle.
        cmap: Colormap for coloring each rectangle.
        norm: Normalization for the colormap.
        rect_height (float): The height (thickness) of each rectangle along the radius.
        zorder (int, optional): The drawing order. Defaults to 50.
    """
        n_segments = len(misc_temp_list_sec)
        angles = np.linspace(angle_start, angle_end, n_segments + 1)
        # rect_height =
        for idx in range(n_segments):
            # Midpoint angle for each rectangle
            angle = (angles[idx] + angles[idx + 1]) / 2
            theta_diff = angles[idx + 1] - angles[idx]

            color = cmap(norm(misc_temp_list_sec[idx]))

            rect = Rectangle(
                xy = (angle - theta_diff, radius - rect_height),  # (angle_start, radius_start)
                width = theta_diff,
                height= rect_height,
                facecolor=color,
                edgecolor='black',
                zorder=0,
                linewidth=0.5
            )

            # Add the rectangle to the plot
            ax.add_patch(rect)

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
            sm, ax=ax, aspect=39, orientation="vertical", pad = 0.1
        )  # Unified colorbar
        if self.type_flag == "misc_T":
            cbar.set_label("$T_{melt}$ - $T_{misc}$ (K)", fontsize=self.fontsize)
        elif self.type_flag == 'e_hull':
            cbar.set_label("$E_{hull}$ (eV/atom)", fontsize=self.fontsize)
        elif self.type_flag == 'melt':
            cbar.set_label("$T_{melt}$ (K)", fontsize=self.fontsize)
        elif self.type_flag == 'entropy':
            cbar.set_label("$-Entropy*k_b$", fontsize=self.fontsize)
        elif self.type_flag == 'elastic':
            cbar.set_label("$E (GPa)$", fontsize=self.fontsize)
        elif self.type_flag == 'density':
            cbar.set_label("$\\rho (g/cc)$", fontsize=self.fontsize)
        elif self.type_flag == 'gibbs':
            cbar.set_label("$G_{mix} (meV/atom)$", fontsize=self.fontsize)

        plt.subplots_adjust(bottom=0.15, top = 0.8)

    def find_misc_temperatures(
        self, member_pos: List[int], N: int, flag: Literal["add", "transmutate"]
    ) -> np.ndarray:
        misc_temp_list = self.misc_temp(member_pos=member_pos, x=self.x, flag=flag, N=N)

        misc_temp_list = np.array([-3000 if np.isnan(i) else i for i in misc_temp_list])
        print(misc_temp_list)
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
        line_colors = ["#009988", "#EE7733", "#0077BB", "#CC3311", "#EE3377", "#BBBBBB"]
        temperature = kwargs.get("temperature", None)

        if self.type_flag == "entropy":
            misc_temp_list = self.entropy(
                member_pos=member_pos,
                x=self.x,
                N=N,
                flag="add",
            )

        if self.type_flag == "elastic":
            misc_temp_list = self.elastic(
                member_pos=member_pos,
                x=self.x,
                N=N,
                flag="add",
            )

        if self.type_flag == "melt":
            misc_temp_list = self.melt_T(
                member_pos=member_pos,
                x=self.x,
                N=N,
                flag="add",
            )

        if self.type_flag == "density":
            misc_temp_list = self.density(
                member_pos=member_pos,
                x=self.x,
                N=N,
                flag="add",
            )

        if self.type_flag == "gibbs":
            misc_temp_list = self.gibbs(
                member_pos=member_pos,
                x=self.x,
                N=N,
                flag="add",
            )

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
            with open("/Users/mcube/Desktop/Projects/far_heaa/src/far_heaa/wenhao/InverseHullWeb-main/InverseHullWeb/tutorial/mol_dict_done.pkl",
                      "rb") as resultFile:
                result = pickle.load(resultFile)

            # misc_temp_list = [result [str]]
            misc_temp_list = [result[str(mol)] for mol in misc_temp_list]
            print(misc_temp_list)
        elif self.type_flag == "misc_T":
            misc_temp_list = self.find_misc_temperatures(member_pos, N, flag="add")


        # self.plot_line(
        #     angle_degrees=angle,
        #     x_values=self.x * pm.distance_calculator(n_alloy, N),
        #     temp_list=misc_temp_list,
        #     ax=ax,
        #     cmap=self.cmap,
        #     norm=self.norm,
        #     zorder=1,
        # )
        self.plot_rectangle(
            angle_degrees=angle,
            x_values=self.x * pm.distance_calculator(n_alloy, N),
            temp_list=misc_temp_list,
            ax=ax,
            cmap=self.cmap,
            norm=self.norm,
            zorder=1,
            bar_width = kwargs.get("bar_width", None)
        )
        angle_radians = np.radians(angle)
        ax.vlines(
            angle_radians,
            ymin=self.y_bias,
            ymax=(self.x * pm.distance_calculator(n_alloy, 1))[-1]
            + self.y_bias,
            linestyles="-",
            color=line_colors[idx2],
            zorder=0,
            alpha=0.7,
            linewidth=1.75,
        )
        rotation = self.text_flipper(angle=angle)
        if '-' not in i:
            ax.text(
                angle_radians,
                pm.distance_calculator(n_alloy, 1) +
            N * 0.19
            + self.y_bias,
                i,
                ha="center",
                va="center",
                color="#882255",
                rotation=rotation,
                weight = 'bold',
                fontsize = self.fontsize
            )
        else:
            if rotation > 360:
                ax.text(
                    angle_radians,
                    pm.distance_calculator(n_alloy, 1)
                +N * 0.20
                # - 0.02 * n_alloy
                + self.y_bias,
                    i,
                    ha="center",
                    va="center",
                    color="black",
                    rotation=rotation,
                    fontsize= self.fontsize
                )
            else:
                ax.text(
                    angle_radians,
                    pm.distance_calculator(n_alloy, 1)
                +N * 0.20
                # - 0.02 * n_alloy
                + self.y_bias,
                    i,
                    ha="center",
                    va="center",
                    color="black",
                    rotation=rotation,
                    fontsize=self.fontsize
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

        if self.type_flag == "e_hull" and temperature is None:
            raise ValueError("Temperature should be provided for e_hull plot.")
        total = []
        for i in MultinaryCombinations.create_multinary(element_list=self.composition, no_comb=list(range(1, len(self.composition)))).values():
            total += list(i)

        total = [i for i in total if self.composition[0] in i]
        comp_str = '-'.join(sorted(self.composition))
        # print(comp_str, total)
        total_mirror = ['-'.join(list(set(self.composition).difference(set(i.split('-'))))) for i in total]

        total += total_mirror
        print(total)
        for idx2, comp in enumerate(total):
            angle = angles[count]
            temp_i = comp.split("-")
            member_pos = fLE.find_indices(self.composition, temp_i)
            if temperature:
                scatter = self.make_one_bar(
                    ax,
                    n_alloy,
                    member_pos,
                    comp,
                    len(temp_i),
                    float(angle),
                    len(temp_i)-1,
                    temperature=temperature,
                    bar_width=len(total)
                )
            else:
                scatter = self.make_one_bar(
                    ax, n_alloy, member_pos, comp, len(temp_i), float(angle), len(temp_i)-1, bar_width=len(total)
                )

            count += 1
        # for idx2, N in enumerate(range(1, len(self.composition))):
        #
        #     combs = (
        #         self.composition
        #         if N == 1
        #         else list(
        #             MultinaryCombinations.create_multinary(
        #                 element_list=self.composition, no_comb=[N]
        #             ).values()
        #         )[0]
        #     )
        #
        #     for idx, i in enumerate(combs):
        #         angle = angles[count]
        #         temp_i = i.split("-")
        #         member_pos = fLE.find_indices(self.composition, temp_i)
        #         if temperature:
        #             scatter = self.make_one_bar(
        #                 ax,
        #                 n_alloy,
        #                 member_pos,
        #                 i,
        #                 N,
        #                 float(angle),
        #                 idx2,
        #                 temperature=temperature,
        #             )
        #         else:
        #             scatter = self.make_one_bar(
        #                 ax, n_alloy, member_pos, i, N, float(angle), idx2
        #             )
        #
        #         count += 1

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
        with open("./mol_dict.pkl", "wb") as resultFile:
            pickle.dump(self.mol_dict, resultFile)
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
            # melt_T.append(self.tm.avg_T_melt(composition=[i], mol_ratio=mol_ratio))

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
        total = pm().total_num_bars(n_alloy)
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
                        bar_width=total+5
                    )
                else:
                    scatter = self.make_one_bar(
                        ax, n_alloy, member_pos, i, N, float(angle), idx2, bar_width=total+5
                    )

                angle_radians = np.radians(angle)
                if transmute_indices:
                    if idx2 != 0 and idx == min(transmute_indices):
                        misc_temp_list_sec = self.gibbs(
                                                member_pos=member_pos,
                                                x=self.x,
                                                N=N,
                                                flag="transmutate",
                                            )
                        self.plot_colored_secant_rectangle(
                            ax = ax,
                            radius=float(
                                (self.x * pm.distance_calculator(n_alloy, n_alloy - 1))[
                                    -1
                                ]
                                + self.y_bias
                            ),
                            angle_start=angle_radians,
                            angle_end=np.radians(angles[transmute_indices[1]] + idx2 * 180),
                            misc_temp_list_sec=misc_temp_list_sec,
                            cmap=self.cmap,
                            norm=self.norm,
                            rect_height=(self.x * pm.distance_calculator(n_alloy, n_alloy - 1))[-1]
                                 - (self.x * pm.distance_calculator(n_alloy, n_alloy - 1))[-2]
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
