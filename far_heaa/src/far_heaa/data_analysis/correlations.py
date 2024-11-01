import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from far_heaa.io.json_handler import JSONHandler

df = pd.read_csv('total_total_features_BCC.csv')
data = np.array(df.values[:, 1:]).astype(float)
cols = list(df.columns[1:])
print(cols)
#normalize data between 0 and 1
# data = (data - np.min(data))/(np.max(data) - np.min(data))


binary_dict = JSONHandler.load_json(folder_path='../database/', file_name='bokas_omegas_processed')
values = list(binary_dict.values())
real_values = np.array([value['BCC'] for value in values])
#normalize values between 0 and 1
# values = (real_values - np.min(real_values))/(np.max(real_values) - np.min(real_values))


sns.histplot(real_values, color='blue', kde=True, label='True values', bins = 50)
plt.xlabel('$\Omega$')
# sns.histplot(data[:, 0], color='red', kde=True, label='Feature1', bins = 50)
# sns.histplot(data[:, 1], color='green', kde=True, label='Feature2', bins = 50)


# plt.scatter(data[:, 0],data[:, 1], c =real_values, cmap ='coolwarm', s=30)
# plt.scatter(data[:,0], real_values)
# plt.scatter(abs(data[:,-1]-data[:,-2]), real_values)
plt.show()