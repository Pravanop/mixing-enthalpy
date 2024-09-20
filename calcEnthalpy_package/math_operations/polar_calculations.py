import numpy as np


class PolarMaths:
	
	@staticmethod
	def polar_to_cartesian(radius, angle):
		x = radius * np.cos(angle)
		y = radius * np.sin(angle)
		return x, y
	
	@staticmethod
	def cartesian_to_polar(x, y):
		r = np.sqrt(x ** 2 + y ** 2)  # Calculate the radius
		theta = np.arctan2(y, x)  # Calculate the angle in radians
		return r, theta
	
	@staticmethod
	def divide_circle_degrees(n):
		# Generate n equally spaced angles in degrees
		return np.linspace(0, 360, n, endpoint=False)
