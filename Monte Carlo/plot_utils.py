import os
import pickle
from collections import Counter
from datetime import datetime
import plotly
import pandas as pd

import plotly.express as px

import matplotlib.pyplot as plt
import numpy as np

def array_plotter(arr, temp, ele_assign):
    z, x, y = arr.nonzero()

    # Extract non-zero elements
    ele = arr[x, y, z].astype(str)

    ele_assign = ele_assign
    inv_eleassign = {str(v): k for k, v in ele_assign.items()}

    color_scheme = {
        'Cr': "#1f77b4",  # Dark Blue
        'W': "#ff7f0e",  # Dark Orange
        'V': "#2ca02c",  # Dark Green
        'Ta': "#d62728",  # Dark Red
        'Ti': "#9467bd"  # Dark Purple
    }
    ele_name = [inv_eleassign[i] for i in ele]

    df = pd.DataFrame({
        'x': x,
        'y': y,
        'z': z,
        'ele': ele_name
    })

    fig = px.scatter_3d(df,
                        x='x',
                        y='y',
                        z='z',
                        color='ele',
                        color_discrete_map=color_scheme,
                        title=f"{temp}",
                        )


    return fig


def retrieve_latest_dump(folder_path="./dump", filter="300"):
    lfolder = os.listdir(folder_path)
    lfolder = [x for x in lfolder if filter in x]
    ans = [datetime.fromisoformat(os.path.splitext(i.split('_')[-1])[0]) for i in lfolder]
    latest = max(ans)
    result = [v for v in lfolder if str(latest) in v]
    print(f"Opening {result[0]}")
    with open(f'./dump/{result[0]}', 'rb') as handle:
        b = pickle.load(handle)

    return b, result[0]


def sro(mol_frac_dict, key, prob):
    """
                                    alpha_ij = 1 - Pij/cj
    """
    return np.round(1 - prob / mol_frac_dict[key[0]], 3)


def sro_list(arr, neighbour_list, ele_assign):
    keys, values = np.unique(arr, return_counts=True)
    ele_dict = dict(zip(keys,values))
    del ele_dict[0]
    total = sum(list(ele_dict.values()))
    ele_dict = {k:v/total for k, v in ele_dict.items()}
    sro_counter = {}
    inv_eleassign = {str(v): k for k, v in ele_assign.items()}
    for point in neighbour_list:
        neighbours = [tuple(sorted((arr[point[0]], arr[i]))) for i in point[1]]
        counter = dict(Counter(neighbours))
        total_bonds = sum(list(counter.values()))
        counter = {key: sro(ele_dict, key, value / total_bonds) for key, value in counter.items()}
        for key, value in counter.items():
            key_modif = '-'.join([inv_eleassign[str(i)] for i in key])
            if key_modif not in sro_counter:
                sro_counter[key_modif] = value
            else:
                sro_counter[key_modif] += value

    sro_counter = {k:round(v/len(neighbour_list),2) for k, v in sro_counter.items()}
    return sro_counter

def line_plot(plot_conf, x, y):

    # Create a new figure
    # plt.figure(figsize=plot_conf['fig_size'])
    fig, ax = plt.subplots(2, 1, figsize=plot_conf['fig_size'])
    # Plot the line with custom style
    ax[0].plot(x, y, color='royalblue', linewidth=2.5, linestyle='-', marker='o', markersize=7, markerfacecolor='orange',
             markeredgewidth=2, alpha=0.8)
    ax[1].hist(y, bins=30, edgecolor="black", color='royalblue', alpha=0.8)
    # Add titles and labels
    ax[0].set_title(plot_conf['title'], fontsize=plot_conf['fontsize'] + 1, fontweight='bold')
    ax[0].set_xlabel(plot_conf['xlabel'], fontsize=plot_conf['fontsize'])
    ax[0].set_ylabel(plot_conf['ylabel'], fontsize=plot_conf['fontsize'])
    ax[1].set_xlabel(plot_conf['ylabel'], fontsize=plot_conf['fontsize'])

    plt.savefig(plot_conf['file_path'], dpi=300)

    # Show the plot
    plt.show()
