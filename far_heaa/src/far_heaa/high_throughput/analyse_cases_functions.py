from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.math_operations.thermo_calculations import ThermoMaths
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import matplotlib.patches as mpatches


def plot_everything(file, data, cases, case, element_list, system):
    tm = ThermoMaths()
    print(f'\n{case}')
    reduction = []
    addition = []
    temp_reduction = []
    temp_addition = []

    h_alloy_addition = []
    h_alloy_reduction = []

    h_comp_addition = []
    h_comp_reduction = []
    reduce_flag = []
    h_alloy_list = []
    h_comp_list = []
    temp_list = []
    for i in cases[case]:
        keys = i.split('.')
        temp_array = np.array(file[keys[0]][keys[1]])
        composition = keys[0].split('-')
        if keys[1] in composition:
            composition.remove(keys[1])
        # temp_array[temp_array == -1000.0] = np.nan
        alloy_temp = tm.avg_T_melt(composition=composition, mol_ratio=[1/system]*system)
        element_temp = tm.avg_T_melt(composition=composition + [keys[1]], mol_ratio=[1/(system+1)]*(system+1))

        h_alloy = 0
        for j in composition:
            h_alloy += data['-'.join(sorted([j, keys[1]]))]['BCC']
        h_alloy_list.append(h_alloy)
        if h_alloy >= 0:
            h_alloy_addition.append(keys[1])
        else:
            h_alloy_reduction.append(keys[1])

        h_comp = 0
        for j in list(MultinaryCombinations.create_multinary(composition, no_comb=[2]).values())[0]:
            h_comp += data[j]['BCC']
        h_comp_list.append(h_comp)
        if h_comp >= 0:
            h_comp_addition.append(keys[1])
        else:
            h_comp_reduction.append(keys[1])

        # if np.mean(temp_array[1:]) >= temp_array[0]:
        # reduction_dist.append(temp_array[-1] - temp_array[0])
        if (temp_array[-1] >= temp_array[0]) or (temp_array[-1] > 0 and np.isnan(temp_array[0])):
            reduction.append(keys[1])
            reduce_flag.append(1)
        else:
            print(temp_array[-1], temp_array[0])
            addition.append(keys[1])
            reduce_flag.append(0)

        temp_list.append(alloy_temp - element_temp)
        if alloy_temp - element_temp >= 0:
            temp_reduction.append(keys[1])
        else:
            temp_addition.append(keys[1])

    fig, ax = plt.subplots(4, 1, sharex=True)
    count_reduction = {key: Counter(reduction).get(key, 0) for key in element_list}
    count_addition = {key: Counter(addition).get(key, 0) for key in element_list}
    colors_main = ['#999933', '#CC6677']
    sns.barplot(y=np.array(list(count_reduction.values())), x=element_list, color=colors_main[0], linewidth=2,
                edgecolor='black', ax=ax[0])
    sns.barplot(y=-np.array(list(count_addition.values())), x=element_list, color=colors_main[1], linewidth=2,
                edgecolor='black', ax=ax[0])
    colors = colors_main  # Colors matching those in the colormap
    labels = ['Increase', 'Decrease']

    # Create custom patches for the legend
    patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, labels)]

    # Add the legend to the plot
    ax[0].legend(handles=patches, title="Miscibility", bbox_to_anchor = (1.01, 0.8), frameon = False)

    count_reduction_temp = {key: Counter(temp_reduction).get(key, 0) for key in element_list}
    count_addition_temp = {key: Counter(temp_addition).get(key, 0) for key in element_list}

    sns.barplot(y=-np.array(list(count_reduction_temp.values())), x=element_list, color=colors_main[1], linewidth=2,
                edgecolor='black', ax=ax[1])
    sns.barplot(y=np.array(list(count_addition_temp.values())), x=element_list, color=colors_main[0], linewidth=2,
                edgecolor='black', ax=ax[1])

    colors = colors_main  # Colors matching those in the colormap
    labels = ['Increase', 'Decrease']

    # Create custom patches for the legend
    patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, labels)]

    # Add the legend to the plot
    ax[1].legend(handles=patches, title="Melt T", bbox_to_anchor=(1.3, 0.8), frameon=False)

    count_reduction_h_alloy = {key: Counter(h_alloy_reduction).get(key, 0) for key in element_list}
    count_addition_h_alloy = {key: Counter(h_alloy_addition).get(key, 0) for key in element_list}
    sns.barplot(y=np.array(list(count_reduction_h_alloy.values())), x=element_list, color=colors_main[1], linewidth=2,
                edgecolor='black', ax=ax[2])
    sns.barplot(y=-np.array(list(count_addition_h_alloy.values())), x=element_list, color=colors_main[0], linewidth=2,
                edgecolor='black', ax=ax[2])
    colors = colors_main[::-1]  # Colors matching those in the colormap
    labels = ['Negative', 'Positive'][::-1]

    # Create custom patches for the legend
    patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, labels)]

    # Add the legend to the plot
    ax[2].legend(handles=patches, title="$H_a$", bbox_to_anchor=(1.01, 0.8), frameon=False)
    count_reduction_h_comp = {key: Counter(h_comp_reduction).get(key, 0) for key in element_list}
    count_addition_h_comp = {key: Counter(h_comp_addition).get(key, 0) for key in element_list}
    sns.barplot(y=-np.array(list(count_reduction_h_comp.values())), x=element_list, color=colors_main[1], linewidth=2,
                edgecolor='black', ax=ax[3])
    sns.barplot(y=np.array(list(count_addition_h_comp.values())), x=element_list, color=colors_main[0], linewidth=2,
                edgecolor='black', ax=ax[3])

    colors = colors_main[::-1]  # Colors matching those in the colormap
    labels = ['Positive', 'Negative'][::-1]

    # Create custom patches for the legend
    patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, labels)]

    # Add the legend to the plot
    ax[3].legend(handles=patches, title="$H_c$", bbox_to_anchor=(1.01, 0.8), frameon=False)

    ax[1].set_ylabel('# Samples')
    ax[0].set_ylabel('# Samples')
    ax[2].set_ylabel('# Samples')
    ax[3].set_ylabel('# Samples')
    plt.subplots_adjust(hspace=0.15, wspace=0, right = 0.80, top = 0.95, bottom=0.1)
    plt.savefig(f'../plots/high_throughput/{case}_{system}.png', dpi = 150)

    print('# Paths that increase miscibility: ', len(reduction))
    print('# Paths that decrease miscibility: ', len(addition))