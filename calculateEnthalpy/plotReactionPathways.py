import pickle
import random
from random import randint

import matplotlib.pyplot as plt
import numpy as np


class PlotReactionPathways:

    def __init__(self, rP, pathway_energies_temp, pathway_scores):
        self.all_pathways = None
        self.rP = rP
        self.pathway_energies_temp = pathway_energies_temp
        self.pathway_scores = pathway_scores

        self.x = {1: 1, 2: 4, 3: 7, 4: 10}
        #color
        self._init_colordata()

        #plotdata
        self.fontsize = 14
        self.linewidth = 1.5
        self.save = True
        self.plot_save_path = "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/plots/"

    def _init_colordata(self):
        self.color = plt.cm.magma(np.linspace(0.0, 1, len(self.rP['all_pathways'])))
        self.link_color = plt.cm.seismic(np.linspace(0.0, 1, len(self.rP['all_pathways'])))
        self.color_dict = dict(zip(list(self.rP['enthalpy_dict'].keys()), self.color))

    @staticmethod
    def _check_length(key: str):
        return len(key.split('-'))

    def level(self, x, color, y, label):
        xlo = 0.5
        xhi = 0.5
        return plt.plot([x - xlo, x + xhi], [y, y], color=color, linestyle='-', linewidth=self.linewidth-1)

    def level_texts(self, x, y, label):
        xlo = 0.5
        xhi = 0.5
        return plt.text(s=label, x=x-xlo, y=y, fontsize=self.fontsize-8, fontweight='bold')

    def link(self, x1, x2, y1, y2, color, best, zorder):
        xlo = 0.5
        xhi = 0.5
        if best:
            return plt.plot([x1 + xhi, x2 - xlo], [y1, y2], color=color, linestyle='-', linewidth=self.linewidth + 0.5,
                            zorder=zorder)
        else:
            return plt.plot([x1 + xhi, x2 - xlo], [y1, y2], color=color, linestyle='--', linewidth=self.linewidth,
                            alpha=0.5, zorder=zorder)

    def plot_texts(self):
        texts = []
        keys = list(self.rP['misc_temp_dict'].keys())
        values = list(self.rP['misc_temp_dict'].values())
        counter_list = []
        for idx, key in enumerate(keys):

            if key not in counter_list:
                close_to = [keys[i] for i in range(len(keys)) if abs(values[i] - values[idx]) <= 150 and self.x[self._check_length(keys[i])] == self.x[self._check_length(key)]]
                for id, i in enumerate(close_to):
                    if i not in counter_list:
                        idx2 = keys.index(i)
                        texts.append(self.level_texts(x=self.x[self._check_length(key)] + id*0.5,
                                                      y=values[idx2]+10,
                                                      label=keys[idx2]))

                        counter_list.append(i)

        return texts

    def plot_levels(self):
        levels = []

        for key, value in self.rP['misc_temp_dict'].items():
            levels.append(self.level(self.x[self._check_length(key)], y=value, color=self.color_dict[key], label=key))

        return levels

    def plot_links(self):

        links = []
        all_paths = list(np.array(list(self.pathway_energies_temp.keys()))[self.pathway_scores[1]])
        self.all_pathways = all_paths
        for idx_path, path in enumerate(all_paths):
            keys = list(self.pathway_energies_temp[path].keys())
            values = list(self.pathway_energies_temp[path].values())
            for idx, link_value in enumerate(values):
                if idx + 1 == len(keys):
                    continue
                if idx_path == 0:
                    best = True
                else:
                    best = False
                links.append(self.link(
                    x1=self.x[self._check_length(keys[idx])],
                    y1=link_value['misc_T'],
                    x2=self.x[self._check_length(keys[idx + 1])],
                    y2=values[idx + 1]['misc_T'],
                    color=self.link_color[idx_path],
                    best=best,
                    zorder=idx_path))

            # count += 1

        return links

    def plot_reaction_pathways(self):

        levels = self.plot_levels()
        links = self.plot_links()
        texts = self.plot_texts()

        for i in range(len(levels)):
            levels[i]
            texts[i]

        for i in range(len(links)):
            links[i]

        plt.xlabel("Reaction Coordinate")
        plt.ylabel("Temperature (K)")
        plt.title(f"Deposition Pathways for {'-'.join(self.rP['ele_list_main'])}")

        if self.save:
            plt.savefig(self.plot_save_path + f"/{'-'.join(self.rP['ele_list_main'])}.png", dpi=300)


with open("/calculateEnthalpy/data/output_data/aziz_bcc_2/pathway_energies_temp.p",
          'rb') as f:
    rP, pathway_energies_temp, pathway_scores = pickle.load(f)
plots = PlotReactionPathways(rP, pathway_energies_temp, pathway_scores)
plots.plot_reaction_pathways()
for idx, paths in enumerate(plots.all_pathways):
    print(idx, ' ', paths)