import numpy as np
from typing import Tuple


class PolarMaths:
	"""
	A class to perform mathematical operations involving polar and cartesian coordinates,
	and to divide a circle into equal angles.

	Methods:
		polar_to_cartesian(radius, angle):
			Converts polar coordinates to cartesian coordinates.

		cartesian_to_polar(x, y):
			Converts cartesian coordinates to polar coordinates.

		divide_circle_degrees(n):
			Divides a circle into `n` equal angles in degrees.
	"""
	
	@staticmethod
	def polar_to_cartesian(radius: float, angle: float) -> Tuple[float, float]:
		"""
		Converts polar coordinates (radius, angle) to cartesian coordinates (x, y).

		Args:
			radius (float): The radius in polar coordinates.
			angle (float): The angle in radians.

		Returns:
			Tuple[float, float]: The cartesian coordinates (x, y) corresponding to the polar coordinates.

		Example:
			x, y = PolarMaths.polar_to_cartesian(radius=5.0, angle=np.pi/4)
			# Converts the polar coordinates (5.0, π/4) to cartesian coordinates.
		"""
		x = radius * np.cos(angle)
		y = radius * np.sin(angle)
		return x, y
	
	@staticmethod
	def cartesian_to_polar(x: float, y: float) -> Tuple[float, float]:
		"""
		Converts cartesian coordinates (x, y) to polar coordinates (radius, angle).

		Args:
			x (float): The x-coordinate in cartesian coordinates.
			y (float): The y-coordinate in cartesian coordinates.

		Returns:
			Tuple[float, float]: The polar coordinates (radius, angle), where the angle is in radians.

		Example:
			r, theta = PolarMaths.cartesian_to_polar(x=3.0, y=4.0)
			# Converts the cartesian coordinates (3.0, 4.0) to polar coordinates.
		"""
		r = np.sqrt(x ** 2 + y ** 2)  # Calculate the radius
		theta = np.arctan2(y, x)  # Calculate the angle in radians
		return r, theta
	
	@staticmethod
	def divide_circle_degrees(n: int) -> np.ndarray:
		"""
		Divides a circle into `n` equal angles in degrees.

		Args:
			n (int): The number of divisions or angles.

		Returns:
			np.ndarray: An array of `n` equally spaced angles in degrees.

		Example:
			angles = PolarMaths.divide_circle_degrees(n=6)
			# Divides the circle into 6 equal angles (0, 60, 120, 180, 240, 300 degrees).
		"""
		return np.linspace(0, 360, n, endpoint=False)
