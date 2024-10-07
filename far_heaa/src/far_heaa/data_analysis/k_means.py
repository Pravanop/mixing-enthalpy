import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
from far_heaa.io.json_handler import JSONHandler
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
#support vector regression
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

#load data
# data = np.load('./total_total_features.npy')
df = pd.read_csv('total_total_features_BCC.csv')
data = np.array(df.values[:, 1:]).astype(float)
cols = list(df.columns[1:])
#normalize data to have mean 0 and std 1
data = (data - np.min(data))/(np.max(data) - np.min(data))


binary_dict = JSONHandler.load_json(folder_path='../database/', file_name='bokas_omegas_processed')
values = list(binary_dict.values())
real_values = np.array([value['BCC'] for value in values])
#normalize values between 0 and 1
values = (real_values - np.min(real_values))/(np.max(real_values) - np.min(real_values))

data_with_y = np.hstack((data, values.reshape(-1, 1)))

sns.histplot(values, color='blue', kde=True, label='True values', bins = 50)
plt.savefig('./analysis_plots/true_values_BCC.png', dpi=200)

#correlation matrix
correlation_matrix = np.round(np.corrcoef(data_with_y.T),2)
fig, ax = plt.subplots(figsize=(13, 13))
sns.heatmap(correlation_matrix,
			annot=True,
			square=True,
			cmap='coolwarm',
			cbar_kws={'label': 'Correlation Coefficient'},
			xticklabels=list(cols)+['Omega'],
			yticklabels=list(cols)+['Omega'])
plt.title('Correlation Matrix')
plt.savefig('./analysis_plots/correlation_matrix_BCC.png', dpi=200)

#random forest
X = data
y = values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
print(y_train.shape, y_test.shape)
# model = RandomForestRegressor(n_estimators=200,
# 							  random_state=42,
# 							  max_depth=15,
# 							  min_samples_split=2,
# 							  min_samples_leaf=1)
model = SVR(kernel='rbf', C=300, gamma=0.4, epsilon=.05)
model.fit(X_train, y_train)
y_pred = model.predict(X_test)

mse = mean_squared_error(y_test, y_pred)
print('normalized Mean Squared Error:', mse)
#unnormalize pred values
y_pred = y_pred*(np.max(real_values) - np.min(real_values)) + np.min(real_values)
y_test = y_test*(np.max(real_values) - np.min(real_values)) + np.min(real_values)
mse_real = mean_squared_error(y_test, y_pred)
print('Mean Squared Error:', mse_real)
fig, ax = plt.subplots(figsize=(10, 10))
sns.histplot(np.sqrt(y_test**2 - y_pred**2), color='blue', kde=True, label='True values',
			 bins = 50)
plt.savefig('./analysis_plots/mse_distribution.png', dpi=200)


# #Feature importance
# importances = model.feature_importances_
# indices = np.argsort(importances)[::-1]
# plt.figure(figsize=(10, 12))
# plt.title("Feature importances")
# plt.bar(range(X_train.shape[1]), importances[indices], align="center")
# plt.xticks(range(X_train.shape[1]), cols, rotation=90)
# plt.xlim([-1, X_train.shape[1]])
# plt.ylabel('Feature importance')
# plt.savefig('./analysis_plots/feature_importance_BCC.png', dpi=200)


