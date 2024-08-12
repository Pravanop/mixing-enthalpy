import matplotlib.pyplot as plt
import numpy as np

temp1 = 300
temp2 = 400
temp3 = 500

x = [1, 4, 7, 10]

for i in range(len(x)):
    plt.plot([x[i]-xlo, x[i]+xhi], [temp1, temp1])
    plt.plot([x[i]-xlo, x[i]+xhi], [temp2, temp2])
    plt.plot([x[i]-xlo, x[i]+xhi], [temp3, temp3])

plt.show()

def level(x, color, y, label):
    xlo = 0.5
    xhi = 0.5
    return (plt.plot([x-xlo, x+xhi], [y, y], color = color, linestyle='-', linewidth=2),
            plt.text(x, y+0.1*y, s = label, va = 'center', ha = 'center'))

def link(x1, x2, y1, y2, color):
    xlo = 0.5
    xhi = 0.5
    return plt.plot([x1+xhi, x2-xlo], [y1,y2], color=color, linestyle = '--', linewidth=2)


# for paths, values in rP.pathway_energies_temp.items():
#     labels = list(values.keys())
#     temperatures = [i['misc_T'] for i in list(values.values())]
#     enthalpies = [i['enthalpy'] for i in list(values.values())]
#     enthalpy_score = rP.enthalpy_scores[count]
#     misc_T_score = rP.misc_T_scores[count]
#     x = []
#     sorted_labels = []
#     for i in labels:
#         sorted_label = '-'.join(sorted(i.split('-')))
#         sorted_labels.append(sorted_label)
#         x.append(rP.coords_dict[sorted_label])
#
#     ax[0].plot(x, enthalpies, label=paths, linewidth=3, linestyle='-', color=color[enthalpy_score], alpha = 0.8, zorder = 0)
#     ax[0].scatter(x, enthalpies, label="_nolegend_", marker = 'o', color="black",s = 100, edgecolors='black', zorder = 1)
#     ax[1].scatter(x, temperatures, label="_nolegend_", marker = 'o', color="black",s = 100, edgecolors='black', zorder = 1)
#     ax[1].plot(x, temperatures, label=paths, color=color[misc_T_score], linewidth=3, linestyle='-', zorder = 0)
#
#     # counter_dict = {}
#     # for i, txt in enumerate(labels):
#     #     sorted_txt = '-'.join(sorted(txt.split('-')))
#     #     if sorted_txt not in counter_dict:
#     #         ax[0].text(s=sorted_txt, x=x[i], y=enthalpies[i]-0.015, fontsize=14, horizontalalignment='center', verticalalignment='bottom')
#     #         ax[1].text(s=sorted_txt, x=x[i], y=temperatures[i]-250, fontsize=14, horizontalalignment='center', verticalalignment='bottom')
#     #         counter_dict[sorted_txt] = 1
#     #     else:
#     #         counter_dict[sorted_txt] += 1
#     count += 1
# fontsize = 14
# plt.subplots_adjust(hspace=0.1, wspace=0.1)
# ax[0].set_ylabel("Mixing Enthalpy (ev/atom)", fontsize = fontsize)
# ax[1].set_ylabel("Miscible Temperature (K)", fontsize = fontsize)
# ax[0].tick_params(axis='both', which='major', labelsize=fontsize)
# ax[1].tick_params(axis='both', which='major', labelsize=fontsize)
# ax[1].set_xticks(ticks = list(rP.coords_dict.values()), labels =list(rP.coords_dict.keys()), rotation=90)
# # ax[1].set_xticklabels(list(rP.coords_dict.values()), rotation=90)
# ax[0].legend(loc='upper right', frameon=False, bbox_to_anchor=(1.15, 1), fontsize='small')
# plt.savefig(f'./plots/reaction_pathways_{"-".join(rP.ele_list_main)}.png', dpi = 300)
#
# # plt.show()