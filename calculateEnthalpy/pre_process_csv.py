"""
I feel it is better that we have the binary pairs as json file. Seems a waste to keep .csv files, and it creates so
many issues while picking pairs. This is just a stopgap to convert other existing datasheets, new ones will be
created in the json format directly
"""

import pandas as pd
from itertools import combinations
import numpy as np
import json

df = pd.read_csv("data/input_data/bcc_bokas.csv")
ele_list = df.columns[1:]
df.set_index('ele', inplace=True)
binaries = list(combinations(ele_list , 2))

results_dict = {}
s = '-'
for idx, binary in enumerate(binaries):
	results_dict[s.join(sorted(binary))] = df[binary[0]][binary[1]] #sort within each pair

results_dict = dict(sorted(results_dict.items(), key=lambda item: item[0])) #sort the whole dictionary

with open('data/input_data/bcc_bokas.json' , 'w') as f:
    json.dump(results_dict, f, ensure_ascii=False, indent=4)

