import numpy as np
import pandas as pd

ele = {
    'Cr': -4.10,
    'W':-8.90,
    'V':-5.31,
    'Ta':-8.10,
    'Ti':-4.85,
    'Hf':-6.14,
}
bin = {
    'Cr-Hf':246,
    'Cr-Ta':108,
    'Cr-Ti':111,
    'Cr-V':-53,
    'Cr-W':108,
    'Hf-Ta':102,
    'Hf-Ti':50,
    'Hf-V':168,
    'Hf-W':122,
    'Ta-Ti':63,
    'Ta-V':64,
    'Ta-W':-61,
    'Ti-V':68,
    'Ti-W':-8,
    'V-W':-67
}

ab = {}
for key, value in bin.items():
    ele_list = key.split('-')
    like_bonds = [ele[i] for i in ele_list]

    # ab[key] = np.round((sum(like_bonds)+value/1000)/2,3)
    ab[key] = value
    for i in ele_list:
        ab[f"{i}-{i}"] = ele[i]

# print(ab)

import json

with open('../data/single_energy.json', 'r') as f:
    d = json.load(f)

print(d)

df = pd.read_csv('../data/binary_BCC.csv')
point5 = np.round(np.array(df['0.5']/24),3)

index = df['Unnamed: 0']

lookup_2 = dict(zip(index, point5))


for key,item in d.items():
    lookup_2[f'{key}-{key}'] = item['dft-energy']

print(lookup_2)
