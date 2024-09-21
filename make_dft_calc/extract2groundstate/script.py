import math

from make_dft_calc.prepareVASPRuns.file_utils import load_json_to_dict

omegas_dict = load_json_to_dict("omegas.json")
omegas = omegas_dict["omegas"]
BCC = omegas["BCC"]
FCC = omegas["FCC"]
HCP = omegas["HCP"]
elements = omegas_dict["elements"]
groundstate_dict = {"As": "BCC",
                    "Sb": "BCC",
                    "Bi": "BCC",
                    "Fe": "BCC",
                    "Cr": "BCC",
                    "W": "BCC",
                    "Mo": "BCC",
                    "V": "BCC",
                    "Ta": "BCC",
                    "Nb": "BCC",
                    "Pb": "FCC",
                    "Ga": "FCC",
                    "Al": "FCC",
                    "In": "FCC",
                    "Cu": "FCC",
                    "Ag": "FCC",
                    "Au": "FCC",
                    "Pd": "FCC",
                    "Pt": "FCC",
                    "Ni": "FCC",
                    "Rh": "FCC",
                    "Ir": "FCC",
                    "Y": "HCP",
                    "Sc": "HCP",
                    "Zr": "HCP",
                    "Hf": "HCP",
                    "Ti": "HCP",
                    "Re": "HCP",
                    "Mn": "HCP",
                    "Os": "HCP",
                    "Ru": "HCP",
                    "Co": "HCP",
                    "Mg": "HCP",
                    "Hg": "HCP",
                    "Cd": "HCP",
                    "Zn": "HCP",
                    "Sn": "HCP",
                    "Ge": "HCP",
                    "Si": "HCP",
                    "Te": "HCP",
                    }
mixing_enthalpies = {}
error_list = []
squared_error = []
x = {}
y = {}
plot_data = {}
for phase in ["BCC", "FCC", "HCP"]:
    mixing_enthalpies[phase] = omegas[phase]
    temp_list = []
    for pair in mixing_enthalpies[phase]:
        a, b = pair.split("-")
        error = 0.5 * (elements[groundstate_dict[a]][a] + elements[groundstate_dict[b]][b]
                                                  - elements[phase][a] - elements[phase][b])
        mixing_enthalpies[phase][pair] = mixing_enthalpies[phase][pair] / 4 \
                                         - 0.5 * (elements[groundstate_dict[a]][a] + elements[groundstate_dict[b]][b]
                                                  - elements[phase][a] - elements[phase][b])
        error_list.append(abs(error))
        temp_list.append(error)
        squared_error.append(math.pow(error, 2))
    plot_data[phase] = temp_list

print(omegas)
import matplotlib.pyplot as plt

# Sample data
index = [i for i in range(1, len(plot_data["BCC"]) + 1)]
error1 = plot_data["BCC"]  # Sample error data set 1
error2 = plot_data["FCC"]  # Sample error data set 2
error3 = plot_data["HCP"]  # Sample error data set 3

# Plotting
plt.figure(figsize=(8, 6))

# Plotting the scatter plots with different colors
plt.scatter(index, error1, label='BCC', color='blue', marker='o')
plt.scatter(index, error2, label='FCC', color='red', marker='s')
plt.scatter(index, error3, label='HCP', color='green', marker='^')

# Customizing the plot
plt.xlabel('Index')
plt.ylabel('Amount of Error')
plt.title('Scatter Plot of Error Data')
plt.legend()
plt.grid(True)

# Show plot
plt.show()