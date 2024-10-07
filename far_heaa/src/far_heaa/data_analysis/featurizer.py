import pandas as pd

from feature_extraction import *
import numpy as np
from far_heaa.io.json_handler import JSONHandler

# elements = ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I']
binary_dict = JSONHandler.load_json(folder_path='../database/', file_name='bokas_omegas_processed')
end_member_data = JSONHandler.load_json(folder_path='../database/', file_name='bokas_end_members_dict')
lattice = 'BCC'
total_total_features = np.zeros((len(list(binary_dict.keys())), 14))
idx = 0
for key, value in binary_dict.items():
	feature_list = []
	elements = key.split('-')
	features, feature_desc = featurizer(elements)
	features_diff = features[:, :7]
	features_diff = abs(features_diff[0] - features_diff[1])
	features_mean = np.mean(features[:,7:], axis=0)
	element_data = [end_member_data[element][lattice] for element in elements]
	feature_list = np.append(features_diff, features_mean)
	feature_list = np.append(feature_list, element_data)
	feature_desc = feature_desc + ['Unary1', 'Unary2']
	
	total_total_features[idx] = feature_list
	idx += 1

df = pd.DataFrame(total_total_features, columns=feature_desc, index = list(binary_dict.keys()))
print(df.head())
df.to_csv(f'total_total_features_{lattice}.csv')
# np.save('./total_total_features', total_total_features)

