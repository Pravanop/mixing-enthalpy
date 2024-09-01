import json
from prepareVASPRuns.file_utils import load_json_to_dict

zhaohan_dict = load_json_to_dict('../data/input_data/bcc_zhaohan.json')
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
                    "Si": "FCC",
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
        mixing_enthalpies[phase][pair] = mixing_enthalpies[phase][pair] / 4 \
                                         - 0.5 * (elements[groundstate_dict[a]][a] + elements[groundstate_dict[b]][b]
                                                  - elements[phase][a] - elements[phase][b])
with open('../data/input_data/bokas2/bcc_bokas.json', 'w') as json_file:
    json.dump(mixing_enthalpies['BCC'], json_file)
with open('../data/input_data/bokasCorrected/fcc_bokas.json', 'w') as json_file:
    json.dump(mixing_enthalpies['FCC'], json_file)
with open('../data/input_data/bokasCorrected/hcp_bokas.json', 'w') as json_file:
    json.dump(mixing_enthalpies['HCP'], json_file)