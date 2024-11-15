# SISSO-based model for predicting DBTT. Needs atomic data from excel spreadsheet. Predicts DBTT in K
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap

#YS
df = pd.read_excel('data/1_Properties_in_Cantor_Space.xlsx')
# dfc = pd.read_excel('data/1_Cantor_benchmark.xlsx')
elements = ['Cr',	'Mn',	'Fe',	'Co',	'Ni']
props = ['Density Avg','FCC YS T C PRIOR','umap0','umap1']
prop = 'entropy'

plt.figure(figsize=(5,4))
composition = df[['Cr', 'Mn', 'Fe', 'Co', 'Ni']].to_numpy()
entropy = -np.sum(composition*np.log(composition + 1e-9), axis = 1)
df['entropy'] = entropy
df = df.sort_values(by=[prop], ascending=False)
plt.scatter(df['umap0'],df['umap1'],c='grey')
# dfp = df[df[prop] > dfc[prop].values[0]]
dfp = df


# cmap = plt.get_cmap('cubehelix', 8)
# cmap = ListedColormap(cmap(np.linspace(0.2, 1, 7))[:-1])
# plt.scatter(dfp['umap0'],dfp['umap1'],c=dfp[prop],cmap=cmap, edgecolor='black')
# plt.colorbar()
# plt.axis('off')
# plt.savefig('results/1_UMAP_RT_YS.png',dpi=500,transparent=True)
# plt.show()

prop = 'Density Avg'
plt.figure(figsize=(5,4))
df = df.sort_values(by=[prop],ascending=False)
plt.scatter(df['umap0'],df['umap1'],c='grey')
# dfp = df[df[prop] < dfc[prop].values[0]]
dfp = df
plt.scatter(dfp['umap0'],dfp['umap1'],c=dfp[prop],cmap='jet',alpha=.7)
plt.colorbar()
plt.axis('off')
plt.savefig('results/1_UMAP_Density.png',dpi=500,transparent=True)
plt.show()
#
#
#
# dfc = pd.read_excel('data/1_Cantor_benchmark.xlsx')
# dfp = pd.read_csv('data/1_Feasible_Alloys.csv')
#
# plt.figure(figsize=(4,4))
# df = df.sort_values(by=[prop])
# plt.scatter(df['umap0'],df['umap1'],c='grey')
# plt.scatter(dfp['umap0'],dfp['umap1'],marker='*',s=300,edgecolors='k',c='Blue')
# plt.scatter(dfc['umap0'],dfc['umap1'],marker='*',s=300,edgecolors='k',c='darkred')
# plt.axis('off')
# plt.savefig('results/1_UMAP_Union_of_Constraints.png',dpi=500,transparent=True)
# plt.show()
#
# elements = ['Cr',	'Mn',	'Fe',	'Co',	'Ni']
# element_colors = {
#     'Co': '#1f77b4',    # Blue
#     'Cr': '#ff7f0e',   # Orange
#     'Fe': '#2ca02c',    # Green
#     'Mn': '#d62728',   # Red
#     'Ni': '#9467bd',   # Purple
# }
#
# plt.figure(figsize=(10, 2))
# for el in elements:
#     sns.distplot(dfp[el], hist=False, label=el,
#                  kde_kws=dict(linewidth=3),color=element_colors[el])
# plt.legend()
# plt.ylim((0, 15))
# plt.xlim((0, 1))
# plt.savefig('results/1_Chemical_Signature.png',transparent=True,dpi=500)
# plt.show()