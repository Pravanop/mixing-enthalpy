import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

#YS
df = pd.read_excel('data/1_Properties_in_Cantor_Space.xlsx')
elements = ['Co','Cr','Fe','Mn','Ni']
df[elements] = (df[elements]*100).astype(int)

order=  list(map(int, range(0, 100, 5)))

element_colors = {
    'Co': '#1f77b4',    # Blue
    'Cr': '#ff7f0e',   # Orange
    'Fe': '#2ca02c',    # Green
    'Mn': '#d62728',   # Red
    'Ni': '#9467bd',   # Purple
}

prop = 'Density Avg'
for el in elements:
    plt.figure(figsize=(5, 2))
    boxprops = dict(facecolor=element_colors[el], edgecolor='black', )
    whiskerprops = dict(color='black')
    medianprops = dict(color='black' )
    capprops = dict(color='black')
    ax = sns.boxplot(x=df[el], y=df[prop], showfliers=False, boxprops=boxprops, whiskerprops=whiskerprops,
                     medianprops=medianprops, capprops=capprops,order=order)
    plt.title(el)
    plt.savefig('results/1_BW_YS_{}.png'.format(el),dpi=500,transparent=True)

prop = 'Density Avg'
for el in elements:
    plt.figure(figsize=(5, 2))
    boxprops = dict(facecolor=element_colors[el], edgecolor='black', )
    whiskerprops = dict(color='black')
    medianprops = dict(color='black' )
    capprops = dict(color='black')
    ax = sns.boxplot(x=df[el], y=df[prop], showfliers=False, boxprops=boxprops, whiskerprops=whiskerprops,
                     medianprops=medianprops, capprops=capprops,order=order)
    plt.title(el)
    plt.savefig('results/1_BW_Density_{}.png'.format(el),dpi=500,transparent=True)