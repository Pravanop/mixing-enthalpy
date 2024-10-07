import itertools
import numpy as np
import pandas as pd


element_list = ['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf']

combs = list(itertools.combinations(element_list, 2))
combs = ['-'.join(sorted(list(i))) for i in combs]

# comp = [12.5, 25, 50, 75, 87.5]


# df_np = np.zeros((len(combs), len(comp)))
# print(df_np.shape)

# df = pd.DataFrame(df_np, columns=comp, index=combs)

# df.to_csv('./dft_output_data/FCC_sub_regular.csv')

bcc = dict(zip(element_list, [0]*len(element_list)))
fcc = dict(zip(element_list, [0]*len(element_list)))
hcp = dict(zip(element_list, [0]*len(element_list)))

master_dict = {
    'BCC': bcc,
    'FCC': fcc,
    'HCP': hcp
}

import json
with open('./dft_output_data/unaries.json', 'w') as f:
    json.dump(master_dict, f, ensure_ascii=True, indent=4)