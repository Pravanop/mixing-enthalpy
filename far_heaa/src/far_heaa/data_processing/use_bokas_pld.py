import pandas as pd
import numpy as np
from far_heaa.io.json_handler import JSONHandler

data = JSONHandler.load_json("../database", file_name="bokas_omegas_processed")
replacement = JSONHandler.load_json("/Users/mcube/Downloads/", file_name="omegas")['omegas']
# replacement = pd.read_csv("/Users/mcube/Downloads/41524_2021_626_MOESM1_ESM.csv", header=None)

# print(replacement.head())
#
# for row, col in replacement.iterrows():
#     if str(col[0]) in data:
#         data[str(col[0])]['BCC'] = float(col[1])
#
#
# JSONHandler.save_json(data, "../database", file_name="bokas_old_omegas_processed")

for lattice in replacement:
    for key, value in replacement[lattice].items():
        if key in data:
            data[key][lattice] = value
            print(key)

JSONHandler.save_json(data, "../database", file_name="bokas_old_omegas_processed")
