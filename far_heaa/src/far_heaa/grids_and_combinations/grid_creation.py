from typing import List, Iterable

import numpy as np


class CompositionGrid:
    """
    A class for generating mole fraction grids and handling high-symmetry and transmutation
    operations for compositions of chemical elements.

    Methods:
            create_mol_grid(n, grid_size):
                    Creates a grid of mole fractions for an n-component system.

            generate_points(n, linspace):
                    Generates points from the provided linspace for (n-1) components.

            filter_valid_points(points):
                    Filters out points whose component sums exceed 1.

            create_high_sym_mol_grid(change_idx, x, n, N):
                    Creates high-symmetry mole fraction grids based on the provided indices.

            create_mol_grid_transmutation(transmutation_indice, n, x):
                    Creates a mole fraction grid with transmutation between two components.
    """

    @staticmethod
    def create_mol_grid(n: int, grid_size: int) -> np.ndarray:
        """
        Create a grid of mole fractions for an n-component system.

        Args:
                n (int): The number of components in the system.
                grid_size (int): The number of grid points to generate per component.

        Returns:
                np.ndarray: A grid of shape (grid_size, n), where each row represents the mole fractions of
                                        the components and all rows sum to 1.

        Raises:
                ValueError: If `n` or `grid_size` is less than 2.

        Example::
        
                grid = CompositionGrid.create_mol_grid(3, 10)
                # grid is a 2D numpy array of mole fractions for a ternary system.
        """
        if n < 2:
            raise ValueError("n must be greater than 1")
        if grid_size < 2:
            raise ValueError("grid_size must be greater than 1")

        linspace = np.linspace(0, 1, grid_size)
        points = CompositionGrid.generate_points(n, linspace)
        valid_points = CompositionGrid.filter_valid_points(points)

        # Compute the last component to ensure the sum equals 1
        last_column = 1 - np.sum(valid_points, axis=1, keepdims=True)
        out = np.hstack([valid_points, last_column])

        # Add homogeneous point where all elements have equal mole fraction
        homogeneous_point = np.full((1, n), 1 / n)
        out = np.vstack([out, homogeneous_point])

        return np.round(out, 4)

    @staticmethod
    def generate_points(n: int, linspace: np.ndarray) -> np.ndarray:
        """
        Generate points from the provided linspace for (n-1) components.

        Args:
                n (int): The number of components in the system.
                linspace (np.ndarray): A 1D array of grid points for each component.

        Returns:
                np.ndarray: A 2D array of shape (grid_size^(n-1), n-1) containing grid points
                                        for the first (n-1) components.
        """
        mesh = np.array(np.meshgrid(*[linspace] * (n - 1)))
        return mesh.T.reshape(-1, n - 1)

    @staticmethod
    def filter_valid_points(points: np.ndarray) -> np.ndarray:
        """
        Filter out points where the sum of the components exceeds 1.

        Args:
                points (np.ndarray): A 2D array of points representing the mole fractions of
                                                         the first (n-1) components.

        Returns:
                np.ndarray: A filtered 2D array containing only the points whose sums are less
                                        than or equal to 1.
        """
        return points[np.sum(points, axis=1) <= 1]

    @staticmethod
    def create_high_sym_mol_grid(
        change_idx: List[int], x: Iterable, n: int, N: int
    ) -> np.ndarray:
        """
        Create a high-symmetry mole fraction grid based on the provided indices.

        Args:
                change_idx (List[int]): Indices of the components to modify in the grid.
                x (Iterable): List of mole fractions to use for the grid.
                n (int): The total number of components in the system.
                N (int): The number of components to symmetrically modify.

        Returns:
                np.ndarray: A 2D array where each row represents a high-symmetry configuration
                                        of the components.
        """
        mol_list = []
        for mol in x:
            addition = np.zeros(n)
            addition += [(1 - mol) / n] * n
            for i in change_idx:
                addition[i] = (1 / N - 1 / n) * mol + (1 / n)
            mol_list.append(addition)

        return np.array(mol_list)

    @staticmethod
    def create_mol_grid_transmutation(
        transmutation_indice: List[int], n: int, x: Iterable
    ) -> np.ndarray:
        """
        Create a mole fraction grid with transmutation between two components.

        Args:
                transmutation_indice (List[int]): Indices of the components undergoing transmutation.
                n (int): The total number of components in the system.
                x ( Iterable): List of mole fractions for transmutation.

        Returns:
                np.ndarray: A 2D array where each row represents a configuration with
                                        transmuted mole fractions between two components.
        """
        mols = []
        for i in x:
            subtract = np.zeros(n)
            subtract += 1 / (n - 1)
            subtract[transmutation_indice[0]] += -1 / (n - 1) + i / (n - 1)
            subtract[transmutation_indice[1]] -= i / (n - 1)
            mols.append(subtract)

        return np.array(mols)
