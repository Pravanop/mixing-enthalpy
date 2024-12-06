import pickle
import numpy as np
with open('../database/intermetallic_database_bokas.pickle', 'rb') as handle:
    im_list = pickle.load(handle)

im_names = []

for j in im_list:
    im_names.append('-'.join(sorted([str(i) for i in j.composition.elements])))

im_names = np.array(im_names)

print(im_list[0])