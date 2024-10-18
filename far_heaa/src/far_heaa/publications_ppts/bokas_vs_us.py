import pandas as pd
import numpy as np
from far_heaa.io.json_handler import JSONHandler
import matplotlib.pyplot as plt


df = pd.read_csv('../database/Outputs_BCC_runs.csv')

end_member_data = {

    "Cr": "-9.52005384",
    "Fe": "-8.25682767",
    "Hf": "-9.74533769",
    "Mo": "-10.93282547",
    "Nb": "-10.10322176",
    "Ta": "-11.81975660",
    "Ti": "-7.69671714",
    "V": "-8.96529006",
    "W": "-12.95591411",
    "Zr": "-8.43952352"

}
binary_dict = JSONHandler.load_json(folder_path='../database/', file_name='bokas_omegas_processed')

us = []
bokas = []
for idx, row in df.iterrows():
    h_dft = row[1:].to_numpy().astype(float)[2]
    end_members = row[0].split('-')
    end_members.reverse()
    end_member_energies = np.array([float(end_member_data[end_member]) for end_member in end_members])

    if not np.isnan(h_dft) :
        print(h_dft)
        h_dft = h_dft - 0.5* end_member_energies[0] - 0.5 * end_member_energies[1]
        us.append(h_dft)
        bokas.append(binary_dict[row[0]]['BCC']/4)

bokas = np.array(bokas)
us = np.array(us)
print(len(us))
print(np.mean(abs(bokas-us)))
plt.hist(abs(bokas-us))
# plt.scatter(us, bokas)
# plt.xlim([-0.8, 0.6])
# plt.ylim([-0.8, 0.6])
plt.show()

