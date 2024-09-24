from typing import Union, Literal, List, Tuple
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from far_heaa.grids_and_combinations import grid_creation
from far_heaa.visualization_toolkit.visualizations import Visualizations


class MatrixHeatmap(Visualizations):
    """
    A class to generate heatmaps for the energy hull (E_hull) based on the composition of a given alloy system.

    Methods:
            find_indices(main_list, subset):
                    Finds the indices of the elements in the subset within the main_list.

            get_ehull_matrix():
                    Calculates and returns the mole fraction grid, energy hull values, and temperature grid.

            plot_ehull_matrix():
                    Plots a heatmap of energy hull values against the mole fraction and temperature grids.

    """

    def __init__(
        self,
        composition: List[str],
        add_ele: Union[List[str], str],
        lattice: Literal["FCC", "BCC", "HCP", "min"],
        meta_data: dict,
        save_flag: bool,
        path_type: Literal["add", "transmutate"],
    ):
        """
        Initializes the MatrixHeatmap class with the alloy composition, lattice type, metadata, and flags for saving and path type.

        Args:
                composition (List[str]): The list of elements in the composition.
                add_ele (Union[List[str], str]): Elements to add or transmutate.
                lattice (Literal['FCC', 'BCC', 'HCP', 'min']): The lattice type.
                meta_data (dict): Metadata containing grid size, folder paths, and flags.
                save_flag (bool): Flag to determine if the plots should be saved.
                path_type (Literal['add', 'transmutate']): The type of path to use ('add' or 'transmutate').
        """
        super().__init__(lattice, meta_data)
        self.x = None
        self.composition = composition
        if isinstance(add_ele, str):
            self.add_ele = [add_ele]
        else:
            self.add_ele = add_ele

        self.total_composition = self.composition + self.add_ele
        self.total_composition = list(set(self.total_composition))
        self.type = path_type

        if self.type == "add":
            self.N = len(self.composition)
            self.n = len(self.total_composition)
            self.n_alloy = len(self.composition)
            self.starting_index = self.find_indices(
                self.total_composition, self.composition
            )
        elif self.type == "transmutate":
            assert len(self.add_ele) == 2
            assert self.add_ele[0] in self.composition
            self.N = len(self.add_ele)
            self.n = len(self.total_composition)
            self.n_alloy = len(self.composition)
            self.transmutation_indices = self.find_indices(
                self.total_composition, self.add_ele
            )
            self.end_composition = [
                ele for ele in self.total_composition if ele != self.add_ele[0]
            ]
        else:
            raise ValueError(
                "Invalid path type. Choose between 'add' or 'transmutate'."
            )

        self.mol_grid_size = 10

        self.save_flag = save_flag

    @staticmethod
    def find_indices(main_list: List[str], subset: List[str]) -> List[Union[int, None]]:
        """
        Finds the indices of the elements in the subset within the main_list.

        Args:
                main_list (List[str]): The main list of elements.
                subset (List[str]): The subset of elements whose indices need to be found.

        Returns:
                List[Union[int, None]]: A list of indices corresponding to the elements in the subset.
                                                                Returns None if an element from the subset is not found in the main_list.
        """
        indices = []
        for value in subset:
            try:
                index = main_list.index(value)
                indices.append(index)
            except ValueError:
                indices.append(None)
        return indices

    def get_ehull_matrix(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Calculates and returns the mole fraction grid, energy hull values, and temperature grid.

        Returns:
                Tuple[np.ndarray, np.ndarray, np.ndarray]: The mole fraction grid, energy hull values, and temperature grid.
        """
        self.x = list(np.linspace(0, 1, self.mol_grid_size))

        mol_grid = []
        if self.type == "add":
            mol_grid = grid_creation.CompositionGrid.create_high_sym_mol_grid(
                change_idx=self.starting_index, x=self.x, n=self.n, N=self.N
            )
        elif self.type == "transmutate":
            mol_grid = grid_creation.CompositionGrid.create_mol_grid_transmutation(
                x=self.x, n=self.n, transmutation_indice=self.transmutation_indices
            )
        mol_grid = mol_grid[::-1]
        mol_grid, e_hulls, temp_grid = self.grid_iterator.e_hull_across_grid(
            composition=self.total_composition,
            lattice=self.lattice,
            mol_grid_size=mol_grid,
            single_temp_flag=False,
            temp_gradation=400,
        )
        return mol_grid, e_hulls, temp_grid

    def plot_ehull_matrix(self) -> Tuple[plt.Axes, plt.Figure]:
        """
        Plots a heatmap of energy hull values against the mole fraction and temperature grids.

        Returns:
                Tuple[plt.Axes, plt.Figure]: The matplotlib axes and figure objects for the plot.
        """
        mol_grid, e_hulls, temp_grid = self.get_ehull_matrix()

        df = pd.DataFrame(e_hulls)
        df.columns = temp_grid
        df = df.T
        df = df.apply(lambda x: x * 1000)  # Convert E_hull values to meV/atom
        positions = []
        for col in df.columns:
            try:
                first_zero_index = df[df[col] == 0].index[0]
            except IndexError:
                first_zero_index = df.index[-1]
            positions.append(first_zero_index)

        cmap = sns.cubehelix_palette(
            start=0.5, rot=-0.61, light=0.98, dark=0.35, hue=1, as_cmap=True
        )
        sns.set_theme(rc={"figure.figsize": (6.4, 6.4)})
        sns.set(font_scale=1.2)
        fig, ax = plt.subplots()
        g = sns.heatmap(
            df,
            yticklabels=np.array(temp_grid).astype(int),
            cmap=cmap,
            xticklabels=np.round(self.x, 2),
            cbar_kws={"label": "$E_{hull}$ (meV/atom)", "shrink": 0.8},
            ax=ax,
        )
        g.set_yticklabels(g.get_yticklabels(), rotation=0)
        g.set_xticklabels(g.get_xticklabels(), rotation=0)
        ax.axes.invert_yaxis()
        ax.axhline(y=0, color="k", linewidth=3)
        ax.axhline(y=df.shape[0], color="k", linewidth=3)
        ax.axvline(x=0, color="k", linewidth=3)
        ax.axvline(x=10, color="k", linewidth=3)

        count_prev = np.where(temp_grid == positions[0])[0][0] + 1
        for idx, i in enumerate(positions):
            idx2 = np.where(temp_grid == i)[0][0] + 1
            if idx2 != count_prev:
                ax.plot([idx, idx], [count_prev, idx2], color="black", linestyle="--")
                count_prev = idx2
            ax.plot([idx, idx + 1], [idx2, idx2], color="black", linestyle="--")

        ax.text(s="-".join(self.composition), y=-1, x=-1)
        if self.type == "transmutate":
            ax.text(s="-".join(self.end_composition), y=-1, x=len(mol_grid) - 0.3)
        if self.type == "add":
            ax.text(s="-".join(self.total_composition), y=-1, x=len(mol_grid) - 0.3)
        ax.set_xlabel("X")
        ax.set_ylabel("T (K)")
        plt.subplots_adjust(bottom=0.15, top=0.9, left=0.14, right=0.98)

        if self.save_flag:
            self.save_figure(
                folders=["heatmap_plots", self.type],
                file_name=f"{''.join(self.composition)}_{''.join(self.end_composition)}",
                fig=fig,
            )
        return ax, fig
