import numpy as np


class CompositionGrid:
	"""
	A class to generate N-1 dimensional grids for alloy compositions.

	Methods
	-------
	create_mol_grid(n: int, grid_size: int) -> np.ndarray:
		Generates a grid of compositions with a linear edge density of grid_size.
	generate_points(n: int, linspace: np.ndarray) -> np.ndarray:
		Generates grid points in N-1 dimensions.
	filter_valid_points(points: np.ndarray) -> np.ndarray:
		Filters points where the sum of mole fractions is <= 1.
	"""
	
	@staticmethod
	def create_mol_grid(n: int, grid_size: int) -> np.ndarray:
		"""
		Generate N-1 dimensional grid with a linear edge density of grid_size.

		Parameters
		----------
		n : int
			Number of elements in the alloy system (must be greater than 1).
		grid_size : int
			Linear point density on edges of the grid (must be greater than 1).

		Returns
		-------
		np.ndarray
			A grid of compositions (N points) with valid mole fractions.

		Raises
		------
		ValueError
			If n is less than 2 or grid_size is less than 2.
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
		Generates the grid points in N-1 dimensions.

		Parameters
		----------
		n : int
			Number of elements in the alloy system.
		linspace : np.ndarray
			An array of evenly spaced values between 0 and 1 for grid points.

		Returns
		-------
		np.ndarray
			A reshaped array of grid points in N-1 dimensions.
		"""
		mesh = np.array(np.meshgrid(*[linspace] * (n - 1)))
		return mesh.T.reshape(-1, n - 1)
	
	@staticmethod
	def filter_valid_points(points: np.ndarray) -> np.ndarray:
		"""
		Filters points where the sum of mole fractions is <= 1.

		Parameters
		----------
		points : np.ndarray
			Array of points in N-1 dimensions.

		Returns
		-------
		np.ndarray
			Filtered array of valid points where the sum of mole fractions is <= 1.
		"""
		return points[np.sum(points, axis=1) <= 1]
	
	@staticmethod
	def create_high_sym_mol_grid(
			change_idx: list[int],
			x: list,
			n: int,
			N: int):
		
		mol_list = []
		for mol in x:
			addition = np.zeros(n)
			addition += [(1 - mol) / n] * n
			for i in change_idx:
				addition[i] = (1 / N - 1 / n) * mol + (1 / n)
			mol_list.append(addition)
		
		return np.array(mol_list)
	
	@staticmethod
	def create_mol_grid_transmutation(transmutation_indice: list,
									  n: int,
									  x: list):
		
		mols = []
		for i in x:
			subtract = np.zeros(n)
			subtract += 1 / (n - 1)
			subtract[transmutation_indice[0]] += -1 / (n - 1) + i / (n - 1)
			subtract[transmutation_indice[1]] -= i / (n - 1)
			mols.append(subtract)
		
		return np.array(mols)
