import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from tqdm import tqdm
from matplotlib.patches import Rectangle, Arrow

from calculateEnthalpy.helper_functions.grid_code import create_multinary

def calculate_temp_for_mol(mol, temp_grid, composition, genre, im_list, pD, temp_composition, ind, equi_mol):
    temp_list = []
    for temp in temp_grid:

        mol_ratio = [equi_mol]*len(composition) + [mol]
        mol_ratio[ind] -= mol
        e_above_hull = pD.find_decomp_products(
                composition=temp_composition,
                mol_ratio=mol_ratio,
                temperature=temp,
                lattice=genre,
                batch_tag=True,
                im=im_list
            )[1]
        if np.isclose(e_above_hull, 0, atol=1e-3):
            temp_list.append(0)
        else:
            temp_list.append(e_above_hull)

    return temp_list

def transmute_ele(composition, el_pre, el_post, pD, genre):
    temp_space = 11
    mol_space = 10

    equi_mol = [1/len(composition)]*len(composition)
    mol_grid = np.round(np.linspace(0.0, equi_mol[0], mol_space), 2)

    # mol_grid = np.round(np.linspace(0.0, 1 / (len(composition) + 1), mol_space), 2)
    temp_composition = composition + [el_post]
    ind = temp_composition.index(el_pre)

    temp_grid = np.linspace(200, 3200, temp_space)
    normalized_mol = np.round((mol_grid - np.min(mol_grid)) / (np.max(mol_grid) - np.min(mol_grid)), 2)
    n_alloy = len(temp_composition)
    all_combs = create_multinary(element_list=temp_composition, no_comb=list(range(2, n_alloy + 1)))
    im_list = []
    for dimensionality, alloy_list in all_combs.items():
        if pD.im_flag:
            im_list += pD.get_intermetallic(alloy_list)

    misc_temp = {}
    for idx, mol in enumerate(tqdm(mol_grid)):
        misc_temp[normalized_mol[idx]] = calculate_temp_for_mol(mol, temp_grid, composition, genre, im_list, pD, temp_composition, ind, equi_mol[1])

    df = pd.DataFrame().from_dict(misc_temp)
    df = df.T
    df.columns = temp_grid
    df = df.T
    df = df.apply(lambda x: x * 1000)
    positions = []
    for col in df.columns:
        try:
            first_zero_index = df[df[col] == 0].index[0]
        except IndexError:
            first_zero_index = df.index[-1]

        positions.append(first_zero_index)

    # print(positions)
    # df = df.iloc[::-1]
    cmap = sns.cubehelix_palette(start=.5, rot=-.61, light=.98, dark=.35, hue=1, as_cmap=True)
    sns.set_theme(rc={'figure.figsize': (11.7, 8.27)})
    sns.set(font_scale=1.4, )
    fig, ax = plt.subplots()
    g = sns.heatmap(df, yticklabels=temp_grid.astype(int), cmap=cmap, square=True,
                     cbar_kws={'label': '$E_{hull}$ meV/atom', "shrink": 0.8}, ax=ax)
    g.set_yticklabels(g.get_yticklabels(), rotation=0)
    ax.axes.invert_yaxis()
    ax.axhline(y=0, color='k', linewidth=3)
    ax.axhline(y=df.shape[1] + 1, color='k', linewidth=3)
    ax.axvline(x=0, color='k', linewidth=3)
    ax.axvline(x=10, color='k', linewidth=3)
    count_prev = np.where(temp_grid == positions[0])[0][0] + 1
    for idx, i in enumerate(positions):
        idx2 = np.where(temp_grid == i)[0][0] + 1
        if idx2 != count_prev:
            ax.plot([idx, idx], [count_prev, idx2], color='black', linestyle='--')
            count_prev = idx2
        ax.plot([idx, idx+1],[idx2, idx2], color='black', linestyle='--')

    composition.remove(composition[ind])
    # ax.set_yticks(rotation=0)
    ax.set_xlabel(f"{el_pre}" + "$_{1-x}$" + f"{el_post}" + "$_{x}$" + f"{''.join(composition)}")
    ax.set_ylabel("T (K)")
    plt.subplots_adjust(bottom=0.2, top=0.9)
    return fig
    # plt.savefig(f'heatmap_{"-".join(temp_composition)}.png')