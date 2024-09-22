import numpy as np
from typing import Tuple, Union


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

            distance_calculator(n, N):
                    Calculates the distance between two values `n` and `N`.

            total_num_bars(n):
                    Calculates the total number of bars based on a given value of `n`.

            angle_assigner(length):
                    Assigns angles based on the length of a sequence.
    """

    def __init__(self):
        pass

    @staticmethod
    def polar_to_cartesian(radius: float, angle: float) -> Tuple[float, float]:
        """
        Converts polar coordinates (radius, angle) to cartesian coordinates (x, y).

        Args:
                radius (float): The radius in polar coordinates.
                angle (float): The angle in radians.

        Returns:
                Tuple[float, float]: The cartesian coordinates (x, y) corresponding to the polar coordinates.
        """
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        return x, y

    @staticmethod
    def cartesian_to_polar(
        x: Union[float, np.ndarray], y: Union[float, np.ndarray]
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Converts cartesian coordinates (x, y) to polar coordinates (radius, angle).

        Args:
                x (float or np.ndarray): The x-coordinate in cartesian coordinates or an array of such.
                y (float or np.ndarray): The y-coordinate in cartesian coordinates or an array of such.

        Returns:
                Tuple[np.ndarray, np.ndarray]: The polar coordinates (radius, angle), where the angle is in radians.
        """
        r = np.sqrt(x**2 + y**2)  # Calculate the radius
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
        """
        return np.linspace(0, 360, n, endpoint=False)

    @staticmethod
    def distance_calculator(n: float, N: float) -> float:
        """
        Calculates the distance between two values `n` and `N` using the formula:
        2 * 1 / sqrt(n - N).

        Args:
                n (float): The first value (numerator).
                N (float): The second value (denominator).

        Returns:
                float: The calculated distance.

        Raises:
                ValueError: If `n` is equal to `N`, which would result in division by zero.
        """
        if n == N:
            raise ValueError("N must be smaller than n")
        return 2 * 1 / np.sqrt(n - N)

    @staticmethod
    def total_num_bars(n: int) -> int:
        """
        Calculates the total number of bars using the formula 2^n - 2.

        Args:
                n (int): The input value used to calculate the total number of bars.

        Returns:
                int: The total number of bars.
        """
        return 2**n - 2

    def angle_assigner(self, length: int) -> np.ndarray:
        """
        Assigns angles based on the length of a sequence. If the length is even, the angles
        are halved.

        Args:
                length (int): The length of the sequence to assign angles to.

        Returns:
                np.ndarray: An array of angles in degrees.
        """
        if length % 2 == 0:
            angles = self.divide_circle_degrees(length) / 2
        else:
            angles = self.divide_circle_degrees(length)

        return angles
