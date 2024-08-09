import pickle
import re

import matplotlib.pyplot as plt
import numpy as np


class Node:

    def __init__(self, name, value):
        self.name = name
        self.alloy = name
        self.temperature = value
        self.node_class = len(name.split('-'))
        self.x = self.node_class * 100
        self.node_type = 0
        self.comp = [1 / self.node_class] * self.node_class

    def __str__(self):
        return self.name


class subNode:

    def __init__(self, name, value):
        self.name = name
        self.temperature = value
        self.comp = re.findall(r"\d+\.\d+", name)
        self.alloy = ''.join([i for i in name if not i.isdigit() and i != '.' and i != ','])
        self.node_type = 1
        self.node_class = len(self.comp)
        self.x = (self.node_class - 1 + float(self.comp[-1]) * self.node_class) * 100

    def __str__(self):
        return self.name


class Link:
    def __init__(self, node1, node2, rank, color, best):
        self.node1 = node1
        self.node2 = node2
        self.x1 = node1.x
        self.x2 = node2.x
        self.y1 = node1.temperature
        self.y2 = node2.temperature
        self.rank = rank
        self.color = color
        self.best = best

    def __str__(self):
        return f'{self.node1.alloy} -> {self.node2.alloy}'


class PlotReactionPathways:

    def __init__(self, rP, constraint: str = 'all'):
        self.rP = rP
        #color
        self._init_colordata()
        self.xlo = 15
        self.xhi = 15
        if constraint != 'all':
            self.constraint = [tuple(i.split('-')) for i in constraint]
        else:
            self.constraint = constraint
        #plotdata
        self.fontsize = 14
        self.linewidth = 1.5
        self.save = True
        self.plot_save_path = "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/plots/"

    def _init_colordata(self):
        self.color = plt.cm.magma(np.linspace(0.0, 1, len(self.rP['all_pathways'])))
        self.link_color = plt.cm.seismic(np.linspace(0.0, 1, len(self.rP['all_pathways'])))
        self.color_dict = dict(zip(list(self.rP['enthalpy_dict'].keys()), self.color))

    def level(self, x, color, y, label):
        return plt.plot([x - self.xlo+1, x + self.xhi-1], [y, y], color=color, linestyle='-', linewidth=self.linewidth + 2,
                        label=label, zorder=25)

    def sub_level(self, x, color, y):
        return plt.scatter(x, y, color=color, marker='o', s= 5, zorder=25)

    def level_texts(self, x, y, label):
        return plt.text(s=label, x=x - self.xlo, y=y, fontsize=self.fontsize - 6, fontweight='bold', zorder=25)

    def link(self, x1, x2, y1, y2, color, best, zorder):
        if best:
            return plt.plot([x1, x2], [y1, y2], color=color, linestyle='-',
                            linewidth=self.linewidth + 0.5,
                            zorder=zorder)
        else:
            return plt.plot([x1, x2], [y1, y2], color=color, linestyle='--',
                            linewidth=self.linewidth,
                            alpha=0.5, zorder=zorder)

    def _create_subnodes(self):
        subnodes = []
        for i in self.rP["sub_misc_temp_dict"].keys():
            subnodes.append((subNode(name=i, value=rP["sub_misc_temp_dict"][i])))

        return subnodes

    def _create_nodes(self):
        nodes = []

        for i in self.rP["misc_temp_dict"].keys():
            nodes.append(Node(name=i, value=rP["misc_temp_dict"][i]))

        return nodes
    def _process_pathways(self):
        all_pathways = []
        print(self.constraint)
        if self.constraint == 'all':
            self.constraint = self.rP["all_pathways"]

        for i in self.rP["all_pathways"]:
            if i in self.constraint:
                temp = []
                for j in range(len(i)):
                    temp.append('-'.join(i[:j + 1]))

                all_pathways.append(temp)

        return all_pathways

    def _create_pathways(self, all_pathways, subnodes, nodes):
        Pathways = []
        for i in all_pathways:
            temp = []
            for j in i:
                for k in nodes:
                    sorted_j = '-'.join((sorted(j.split('-'))))
                    if sorted_j == k.alloy:
                        temp.append(k)

                for k in subnodes:
                    sorted_j = '-'.join((sorted(j.split('-'))))
                    if sorted_j == k.alloy:
                        temp.append(k)

            temp = sorted(temp, key=lambda y: y.x)
            Pathways.append(temp)

        return Pathways

    def _rank_pathways(self, pathways):
        pathway_scores = []
        for i in pathways:
            temp_score = 0
            for node in i:
                temp_score += node.temperature

            pathway_scores.append(temp_score/len(i))

        return pathway_scores, np.argsort(pathway_scores)

    def _create_links(self, pathways, order):
        links = []
        for idx in range(len(pathways)):
            best = False
            pathway = pathways[order[idx]]
            if idx == 0:
                best = True
            color = self.link_color[idx]
            for first, second in zip(pathway, pathway[1:]):
                links.append(Link(first, second, idx, color, best))

        return links

    def display_pathways(self, all_pathways, order):
        for ind, idx in enumerate(order):
            print(f"{ind}  {all_pathways[idx]}  {self.scores[idx]}")
    def main_plot(self):

        subnodes = self._create_subnodes()
        nodes = self._create_nodes()
        all_pathways = self._process_pathways()

        pathways = self._create_pathways(all_pathways, subnodes, nodes)
        self.scores, order = self._rank_pathways(pathways)

        links = self._create_links(pathways, order)

        for idx, node in enumerate(nodes):
                self.level(x=node.x,
                           color=self.color_dict[node.alloy],
                           y=node.temperature,
                           label=node.name)
                self.level_texts(node.x, node.temperature+10, node.alloy)

        for subnode in subnodes:
                self.sub_level(x=subnode.x,
                               y=subnode.temperature,
                               color="black")

        for link_o in links:
            self.link(x1=link_o.x1,
                      x2=link_o.x2,
                      y1=link_o.y1,
                      y2=link_o.y2,
                      color=link_o.color,
                      best=link_o.best,
                      zorder=link_o.rank)

        plt.ylim([200, 3000])
        plt.xlabel("Reaction Coordinate")
        plt.ylabel("Temperature (K)")
        plt.title(f"Deposition Pathways for {'-'.join(self.rP['ele_list_main'])}")

        self.display_pathways(self.rP["all_pathways"], order)

        if self.save:
            plt.savefig(self.plot_save_path + f"/{'-'.join(self.rP['ele_list_main'])}.png", dpi=300)



with open(
        "/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/data/output_data/aziz_bcc_2/pathway_energies_temp.p",
        'rb') as f:
    rP, pathway_energies_temp, pathway_scores = pickle.load(f)

# plot_pathways = PlotReactionPathways(rP, constraint=['V-Ti-Zr-Nb', 'V-Ti-Nb-Zr', 'V-Nb-Zr-Ti'])
plot_pathways = PlotReactionPathways(rP, constraint="all")
plot_pathways.main_plot()
