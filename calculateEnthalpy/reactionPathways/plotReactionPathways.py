import pickle
import re
from typing import Union, List

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PathCollection
from matplotlib.lines import Line2D
from matplotlib.text import Text


class Node:
    """
    A main node object for equimolar alloys
    Attributes
    name: alloy composition of the form 'El1-El2'
    alloy: alias of name
    temperature: miscibility temperature
    node_class: dim of node
    x: x coordinate of the node for plotting
    node_type: bool that is 0 for equimolar node and 1 for off-equimolar node
    comp: list of mole fractions for the equimolar alloy
    """

    def __init__(self,
                 name: str,
                 value: float):
        """

        Args:
            name:  alloy composition of the form 'El1-El2' (Strongly recommended to alphabetically sort)
            value: miscibility temperature of the alloy
        """
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
    """
    A sub node object for off-equimolar alloys
    Attributes
    name: alloy composition of the form 'El1,x1-El2,x2'
    alloy: strips the composition and creates a name of form 'El1-El2'
    temperature: miscibility temperature
    node_class: dim of node
    x: x coordinate of the node for plotting
    node_type: bool that is 0 for equimolar node and 1 for off-equimolar node
    comp: list of mole fractions for the equimolar alloy
    """

    def __init__(self,
                 name: str,
                 value: float):
        """

        Args:
            name:  alloy composition of the form 'El1,x1-El2,x2' (Strongly recommended to alphabetically sort)
            value: miscibility temperature of the alloy
        """
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
    """
    A link object to connect nodes belonging in the same pathway
    Attributes
    node1: the start node
    node2: the end node
    x1: the x coordinate of the start node
    x2: the x coordinate of the end node
    y1: the miscibility temperature of the start node
    y2: the miscibility temperature of the end node
    rank: the rank of the link based on the pathway score
    color: color of link which is based on which pathway it is
    best: if the link is part of the lowest energy pathway
    """

    def __init__(self,
                 node1: Union[Node, subNode],
                 node2: Union[Node, subNode],
                 rank: int,
                 color: str,
                 best: bool):
        """

        Args:
            node1: Start node
            node2: End Node
            rank: rank of the link based on pathway score
            color: color for link
            best: if the link is part of the lowest score pathway
        """
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
    """
    A class to plot the reaction pathways using matplotlib

    Attributes
    xlo, xhi: length of node line
    fontsize: font size of plot
    save: bool to decide if figure should be saved
    plot_save_path: folder path to save the plot
    """

    def __init__(self,
                 rP: dict,
                 constraint: str = 'all'):
        """
        The pickled dict computed from reactionPathways class
        Args:
            rP: dictionary containing all the information of the deposition pathways computed from the reactionPathways class
            constraint: if only a subset of paths are needed to be plotted
        """
        self.rP = rP
        self._init_colordata()
        self.xlo = 15
        self.xhi = 15
        if constraint != 'all':
            self.constraint = [tuple(i.split('-')) for i in constraint]
        else:
            self.constraint = constraint
        self.fontsize = 14
        self.linewidth = 1.5
        self.save = True
        self.plot_save_path = "/plots/"

    def _init_colordata(self) -> None:
        """
        Initialize the color data based on color maps for nodes and links.
        Returns: None

        """
        self.color = plt.cm.magma(np.linspace(0.0, 1, len(self.rP['all_pathways'])))
        self.link_color = plt.cm.seismic(np.linspace(0.0, 1, len(self.rP['all_pathways'])))
        self.color_dict = dict(zip(list(self.rP['enthalpy_dict'].keys()), self.color))

    def level(self,
              x: float,
              color: str,
              y: float,
              label: str) -> list[Line2D]:
        """
        Plots one node as a short, straight line
        Args:
            x: the x coordinate of the node for plotting
            color: the assigned color for the node
            y: the miscibility temperature of the node
            label: the name of the node

        Returns: list of Line2D objects depicting one node.

        """
        return plt.plot([x - self.xlo + 1, x + self.xhi - 1], [y, y], color=color, linestyle='-',
                        linewidth=self.linewidth + 2,
                        label=label, zorder=25)

    def sub_level(self,
                  x: float,
                  color: str,
                  y: float) -> PathCollection:
        """
        Plots one node as a circle.
        Args:
            x: the x coordinate of the node for plotting
            color: the assigned color for the node
            y: the miscibility temperature of the node

        Returns: scatter plot object for each subnode

        """
        return plt.scatter(x, y, color=color, marker='o', s=5, zorder=25)

    def level_texts(self,
                    x: float,
                    y: float,
                    label: str) -> Text:
        """
        Plots label for each main node. TODO Still needs to resolve overlapping
        Args:
            x: x coordinate of the node for plotting
            y: y coordinate of the node for plotting
            label: text to be plotted, name of node

        Returns: text object

        """
        return plt.text(s=label, x=x - self.xlo, y=y, fontsize=self.fontsize - 6, fontweight='bold', zorder=25)

    def link(self,
             x1: float,
             x2: float,
             y1: float,
             y2: float,
             color: str,
             best: bool,
             zorder: int) -> list[Line2D]:
        """
        Plots a link between two nodes as a dashed line or straight line for best pathway
        Args:
            x1: x coordinate of the start node
            x2: x coordinate of the end node
            y1: y coordinate of the start node
            y2: y coordinate of the end node
            color: color assigned to pathway
            best: if link is part of the lowest energy pathway
            zorder: rank of the pathway as an int

        Returns: list of Line2D objects depicting one link.

        """
        if best:
            return plt.plot([x1, x2], [y1, y2], color=color, linestyle='-',
                            linewidth=self.linewidth + 0.5,
                            zorder=zorder)
        else:
            return plt.plot([x1, x2], [y1, y2], color=color, linestyle='--',
                            linewidth=self.linewidth,
                            alpha=0.5, zorder=zorder)

    def _create_subnodes(self) -> list[subNode]:
        """
        Creates all subnodes
        Returns: list of subnodes

        """
        subnodes = []
        for i in self.rP["sub_misc_temp_dict"].keys():
            subnodes.append((subNode(name=i, value=rP["sub_misc_temp_dict"][i])))

        return subnodes

    def _create_nodes(self) -> list[Node]:
        """
        Creates all nodes
        Returns: list of nodes

        """
        nodes = []

        for i in self.rP["misc_temp_dict"].keys():
            nodes.append(Node(name=i, value=rP["misc_temp_dict"][i]))

        return nodes

    def _process_pathways(self) -> list:
        """
        Creates all reaction pathways from permutations computed in previous module
        Returns: a list such as ['E1', 'E1-E2', 'E1-E2-E3'], with each string sorted alphabetically

        """
        all_pathways = []
        print(self.constraint)
        if self.constraint == 'all':
            self.constraint = self.rP["all_pathways"]

        for i in self.rP["all_pathways"]:
            if i in self.constraint: # to plot certain pathways only
                temp = []
                for j in range(len(i)):
                    temp.append('-'.join(i[:j + 1]))

                all_pathways.append(temp)

        return all_pathways

    def _create_pathways(self,
                         all_pathways: list,
                         subnodes: list,
                         nodes: list) -> list:
        """
        Adds the actual nodes and subnodes to the aliases created by _process_pathways
        Args:
            all_pathways:  a list such as ['E1', 'E1-E2', 'E1-E2-E3'], with each string sorted alphabetically
            subnodes: list of all subnodes
            nodes: list of nodes

        Returns: list of pathways with nodes and sub-nodes in increasing composition and dimensionality

        """
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

    def _rank_pathways(self,
                       pathways: list):
        """
        Ranks the pathways based on the provided pathway scores
        Args:
            pathways: list of all pathways

        Returns:

        """
        pathway_scores = []
        for i in pathways:
            temp_score = 0
            for node in i:
                temp_score += node.temperature

            pathway_scores.append(temp_score / len(i))

        return pathway_scores, np.argsort(pathway_scores)

    def _create_links(self,
                      pathways: list,
                      order: list):
        """
        Creates the links between nodes, subnodes in the pathways
        Args:
            pathways: all pathways
            order: list of ranks

        Returns: list of links

        """
        links = []
        for idx in range(len(pathways)):
            best = False
            pathway = pathways[order[idx]]
            if idx == 0: #best rank
                best = True
            color = self.link_color[idx]
            for first, second in zip(pathway, pathway[1:]):
                links.append(Link(first, second, idx, color, best))

        return links

    def display_pathways(self,
                         all_pathways: list,
                         order: list) -> None:
        """
        Prints the pathways with their scores
        Args:
            all_pathways: all pathways
            order: list of ranks

        Returns: None

        """
        for ind, idx in enumerate(order):
            print(f"{ind}  {all_pathways[idx]}  {self.scores[idx]}")

    def main_plot(self) -> None:
        """
        Main plotting function that plots all nodes, subnodes, links and saves the figure in an alloted folder
        Returns:

        """

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
            self.level_texts(node.x, node.temperature + 10, node.alloy)

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

#example usage
if __name__ == "__main__":
    with open(
            "/data/output_data/aziz_bcc_2/pathway_energies_temp.p",
            'rb') as f: #loading the previous computed pickle file.
        rP, pathway_energies_temp, pathway_scores = pickle.load(f)

    # plot_pathways = PlotReactionPathways(rP, constraint=['V-Ti-Zr-Nb', 'V-Ti-Nb-Zr', 'V-Nb-Zr-Ti'])
    plot_pathways = PlotReactionPathways(rP, constraint="all")
    plot_pathways.main_plot()
