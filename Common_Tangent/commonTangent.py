import itertools
from typing import Union
from scipy.interpolate import CubicSpline
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sympy import Symbol, nsolve
import mpmath
import sympy as sp
from scipy.spatial import ConvexHull


class commonTangent:
    """"
    Finds the set of common tangents for a given set of curves and points
    """

    def __init__(self, omega: Union[dict, float] = 0):

        if isinstance(omega, dict):
            """
            dictionary should be of the form:
            {
                "bcc" : omegai,
                "fcc" : omegaj
                "intermetallic": [h_i, enthalpy in eV]
            }
            """
            self.omega_dict = omega
        else:
            self.omega = omega
        self.k_b = 8.617333262e-05  #eV
        self.temperature = [300, 3000, 100]

    def binary_regular_model_enthalpy(
            self,
            x_i: float,
            omega: float
    ) -> float:
        """
        Basic function to calculate a regular model. Taken from John Cavin's thesis Eq 3.2.
        :param omega:
        :param x_i: mol fraction
        :return: mixing enthalpy value for mole fraction
        """
        return omega * x_i * (1 - x_i)

    def binary_subregular_model_enthalpy(
            self,
            x_i: float,
            omega: list[float, float]
    ) -> float:
        """
                Basic function to calculate a sub-regular model. Taken from John Cavin's thesis Eq 3.3.
                :param omega:
                :param x_i: mol fraction
                :return: mixing enthalpy value for mole fraction
                """
        return x_i * (1 - x_i) * (omega[0] * x_i + omega[1] * (1 - x_i))

    def deriv_binary_regular_model_enthalpy(
            self,
            x_i: float,
            omega: float
    ) -> float:
        """
        Derivative for a binary regular model.
        :param omega:
        :param x_i: mol fraction

        :return: derivative mixing enthalpy value for mole fraction
        """
        return (1 - 2 * x_i) * omega

    def deriv_binary_subregular_model_enthalpy(self,
                                               x_i: float,
                                               omega: list[float, float]) -> float:

        return x_i * (1 - x_i) * (omega[0] - omega[1]) + (1 - 2 * x_i) * (omega[0] * x_i + omega[1] * (1 - x_i))

    def configEntropy(self,
                      x_i: float) -> float:
        """
        A simple function to calculate boltzmann configurational entropy.

                        delta_S = -k_b*sum_i=1_n(x_i*ln(x_i))

        :param x_i: mol fraction
        :return: boltzmann configurational entropy
        """
        return -self.k_b * (x_i * sp.log(x_i) + (1 - x_i) * sp.log(1 - x_i))

    def deriv_binaryconfigEntropy(self,
                                  x_i: float) -> float:
        """

        :param mol_ratio:
        :return:
        """
        return -self.k_b * (sp.log(x_i) - (sp.log(1 - x_i)))

    def gibbs(self,
              x_i,
              omega: Union[float, list[float, float]],
              temperature: float,
              flag: int = 1,
              ) -> float:
        """
        Function to implement G = H - TS
        :param enthalpy: the mixing enthalpy
        :param entropy: the configurational entropy
        :param temperature: the desired temperature
        :return: gibbs free energy
        """
        if flag == 1:
            return self.binary_regular_model_enthalpy(x_i=x_i,
                                                  omega=omega) - temperature * self.configEntropy(
            x_i=x_i)
        if flag == 2:
            return self.binary_subregular_model_enthalpy(x_i=x_i,
                                                         omega=omega) - temperature * self.configEntropy(
                x_i=x_i)

    def deriv_gibbs(self,
                    x_i: float,
                    omega: Union[float, list[float, float]],
                    temperature: float,
                    flag: int = 1
                    ) -> float:
        """
        Function to implement G = H - TS
        :param enthalpy: the mixing enthalpy
        :param entropy: the configurational entropy
        :param temperature: the desired temperature
        :return: gibbs free energy
        """
        if flag == 1:
            return self.deriv_binary_regular_model_enthalpy(x_i=x_i,
                                                        omega=omega) - temperature * self.deriv_binaryconfigEntropy(
            x_i=x_i)
        if flag == 2:
            return self.deriv_binary_subregular_model_enthalpy(x_i=x_i,
                                                               omega=omega) - temperature * self.deriv_binaryconfigEntropy(
                x_i=x_i)

    def solver_curves(self, omega1,omega2, temperature, flag=1):
        mpmath.mp.dps = 15
        x1, x2 = Symbol('x1'), Symbol('x2')
        f1 = (self.gibbs(omega=omega1, temperature=temperature, x_i=x2, flag=flag) - self.gibbs(omega=omega2,
                                                                                                temperature=temperature,
                                                                                                x_i=x1, flag=flag)) / (
                         x2 - x1) - self.deriv_gibbs(x_i=x1, omega=omega2, temperature=temperature, flag=flag)
        f2 = self.deriv_gibbs(x_i=x2, omega=omega1, temperature=temperature, flag=flag) - self.deriv_gibbs(x_i=x1,
                                                                                                           omega=omega2,
                                                                                                           temperature=temperature,
                                                                                                           flag=flag)
        try:
            mat = nsolve((f1, f2), (x1, x2), (0.1, 0.9), tol = 1e-9)
            x1, x2 = mat[0], mat[1]
        except ZeroDivisionError as e:
            x1, x2 = 1.1, 1.1

        return x1, x2

    def solver_points(self, omega, point, temperature, flag = 1, initial = 0.1):
        x1 = Symbol('x1', real=True)
        f1 = (self.gibbs(omega=omega, temperature=temperature, x_i=x1, flag=flag) - point[1])/(x1- point[0]) - self.deriv_gibbs(x_i=x1, omega=omega, temperature=temperature, flag=flag)
        try:
            mat = nsolve((f1), (x1), (initial), tol = 1e-9)
            x1, x2 = mat, point[0]
        except ZeroDivisionError as e:
            x1, x2 = 1.1, 1.1

        return x1, x2

    def multi_temp_solver(self, flag=1):
        temperature_range = np.linspace(self.temperature[0], self.temperature[1], self.temperature[2] + 1)

        solution_dict = {}
        for idx, temp in enumerate(temperature_range):

            keys = list(self.omega_dict.keys())
            combinations = list(itertools.combinations_with_replacement(keys, r=2))
            combinations.remove(('im', 'im'))
            solution_dict[temp] = []
            for idx, comb in enumerate(combinations):
                if 'im' in comb:

                    point = self.omega_dict['im']
                    list(comb).remove('im')
                    omega1 = self.omega_dict[comb[0]]
                    omega2 = omega1
                    sol = self.tangent_single_temp(omega1=omega1,omega2=omega2,point=point,flag=flag, temp=temp)
                    for x in sol:
                        if x != 1.1:
                            y1 = self.gibbs(x_i = x, omega=omega1, temperature=temp, flag=flag)
                            solution_dict[temp].append([x,y1])

                    solution_dict[temp].append(point)

                else:
                    point = None
                    omega1 = self.omega_dict[comb[0]]
                    omega2 = self.omega_dict[comb[1]]
                    sol = self.tangent_single_temp(omega1=omega1, omega2=omega2, point=point, flag=flag, temp=temp)
                    for x in sol:
                        if x != 1.1:
                            y = self.gibbs(x_i=x, omega=omega1, temperature=temp, flag=flag)
                            solution_dict[temp].append([x, y])

            if len(solution_dict[temp]) <= 2:
                solution_dict[temp] = [i[0] for i in solution_dict[temp]]

            else:
                hull = ConvexHull(solution_dict[temp])
                solution_dict[temp] = [solution_dict[temp][i][0] for i in hull.vertices]
                # solution_dict[temp].sort()

        return solution_dict

    def tangent_single_temp(self, point = None, flag = 1, omega1 = None, omega2 = None, temp = None, solution_dict = {}):

        sol = []
        if not point:

            x1, x2 = self.solver_curves(omega1=omega1, omega2=omega2, temperature=temp, flag=flag)
            x1, x2 = round(x1, 4), round(x2, 4)
            if x1 != x2 and x1 != 1.1:
                if isinstance(x2, float):
                    if x1.is_real:
                        if 0.0 <= x1 <= 1.0 and 0.0 <= x2 <= 1.0:
                            sol.append(x1)
                            sol.append(x2)
                else:
                    if x1.is_real and x2.is_real:
                        if 0.0 <= x1 <= 1.0 and 0.0 <= x2 <= 1.0:
                            sol.append(x1)
                            sol.append(x2)
            else:
                sol.append(1.1)
                sol.append(1.1)

        else:
            for i in [0.1, 0.9]:
                x1, x2 = self.solver_points(omega=omega1, temperature=temp, flag=flag, point=point, initial=i)

                x1, x2 = round(x1, 4), round(x2, 4)
                if x1 != x2 and x1 != 1.1:
                    if isinstance(x2, float):
                        if x1.is_real:
                            if 0.0 <= x1 <= 1.0 and 0.0 <= x2 <= 1.0:
                                sol.append(x1)
                                sol.append(x2)
                    else:
                        if x1.is_real and x2.is_real:
                            if 0.0 <= x1 <= 1.0 and 0.0 <= x2 <= 1.0:
                                sol.append(x1)
                                sol.append(x2)
                else:
                    sol.append(1.1)
                    sol.append(1.1)

        return sol



    def multi_temp_plotter(self, omega1, omega2):
        temperature_range = np.linspace(self.temperature[0], self.temperature[1], self.temperature[2] + 1)
        x_range = np.linspace(0, 1, 101)
        for idx, temp in enumerate(temperature_range):
            if idx % 10 == 0:
                if omega1 == omega2:
                    gibbs1 = [self.gibbs(i, omega1, temp) for i in x_range]
                    plt.plot(x_range, gibbs1, label=f"{temp} K", linewidth=2)
                else:
                    gibbs1 = [self.gibbs(i, omega1, temp) for i in x_range]
                    gibbs2 = [self.gibbs(i, omega2, temp) for i in x_range]
                    plt.plot(x_range, gibbs1, label=f"{temp} K", linewidth=2)
                    plt.plot(x_range, gibbs2, label=f"{temp} K", linewidth=2)

        plt.legend(loc = 'right', bbox_to_anchor=(1.2, 0.5), borderaxespad=0.)


        plt.show()




cT = commonTangent(omega={
    'bcc': 0.108*4,
    'im': [0.3, -0.02]
})

sol = cT.multi_temp_solver(flag=1)
cT.multi_temp_plotter(omega2=0.108*4, omega1=0.108*4)
print("solving done")

points = []
for key, value in sol.items():
    if len(value) >1:
        for k in value:
            points.append([k, key])
    else:
        points.append([value[0], key])

points = np.array(points)

df = pd.DataFrame({'x': points[:, 0], 'y': points[:, 1]})

df.sort_values(by='x', inplace=True)
df.drop_duplicates(subset=['x'], inplace=True)
print(df.head())

cs = CubicSpline(df['x'], df['y'])
x_range = np.linspace(0, 1, 1001)

plt.plot(x_range, cs(x_range), label='Cubic Spline', linewidth=2)
plt.scatter(df['x'], df['y'], color='blue', alpha=0.3)
plt.show()
