import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


unaries = {
"Cr": -9.52540098,
"Fe": -8.26088546,
"Hf": -9.74042020,
"Mo": -10.93737096,
"Nb": -10.09640236,
"Ta": -11.78372388,
"Ti": -7.69369400,
"V": -8.94439489,
"W": -12.96581509,
"Zr": -8.43331072
}

dft_enthalpies = pd.read_csv('/Users/mcube/Downloads/Outputs_BCC_runs.csv')

x_point = dft_enthalpies.columns.to_numpy()[1:].astype(float)/100

arr = dft_enthalpies.to_numpy()


def sub_regular_model(x_i, omega1, omega2):
    return x_i * (1 - x_i) * (omega1 * x_i + omega2 * (1 - x_i))

x_real = np.append(0, x_point)
x_real = np.append(x_real, 1)
systems = {}
for i in arr:
    system = i[0]
    eles = system.split('-')
    h_dft = i[1:]

    y_real = h_dft - x_point*unaries[eles[1]] - (1-x_point)*unaries[eles[0]]
    y_real = np.pad(y_real, 1)
    y_real *= 1000

    popt, pconv = curve_fit(sub_regular_model, xdata = x_real, ydata=y_real)
    systems[system] = [popt, y_real]


x = np.linspace(0, 1, 100)


for key, value in systems.items():
    reg_y_pred = sub_regular_model(x, *value[0])
    plt.plot(x, reg_y_pred)
    plt.scatter(x_real, value[1])
    plt.title(key)
    plt.show()

