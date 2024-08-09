import datetime
from itertools import combinations
from initial_config import InitialConfig
import numpy as np
import random
from tqdm import tqdm
import pickle
from collections import Counter
from nearest_neighbour import create_neighbor_list
from plot_utils import line_plot
import os
from lookup import enthalpy_model_lookup


class MonteCarlo:

    def __init__(self, initial_conf, config_dict):

        self.initial_config = initial_conf
        self.lookup = self.initial_config.lookup
        self.config_dict = config_dict
        self.log = config_dict['log']
        self.energy_trajectory = []
        self.structure_trajectory = []
        self.steps = []
        self.kb = 8.612e-5
        self.neighbour_list = create_neighbor_list(self.initial_config.final_bcclattice, flag=1)
        self.neighbour_list2 = create_neighbor_list(self.initial_config.final_bcclattice, flag=2)
        self.neighbour_dict = {i[0]: i[1] for i in self.neighbour_list}
        self.neighbour_dict2 = {i[0]: i[1] for i in self.neighbour_list2}
    def hamiltonian_bonds(self, point, arr):
        """
        Bonds
        :param point:
        :param arr:
        :return:
        """
        site = arr[point[0]]
        return sum([self.lookup[str(sorted([site, arr[i]]))] for i in point[1]])

    def hamiltonian_enthalpy(self, point, arr):
        """
        Composition
        :param point:
        :param arr:
        :return:
        """
        counter = dict(Counter([arr[i] for i in point[1]] + [arr[point[0]]]))
        total = sum(counter.values())
        counter = {key: round(value / total, 2) for key, value in counter.items()}
        combs = list(combinations(counter.keys(), 2))
        return sum([self.lookup[str(sorted([i[0], i[1]]))] * counter[i[0]] * counter[i[1]] * 4 for i in combs])

    def energy_finder_new(self, arr, neighbour_list, flag):
        if flag == "bonds":
            return sum([self.hamiltonian_bonds(point, arr) for point in neighbour_list]) / 2/2000
        if flag == "enthalpy":
            return sum([self.hamiltonian_enthalpy(point, arr) for point in neighbour_list]) / 8/2000

    def pair_swapper(self, arr, new=False):

        rand = random.choices(self.initial_config.non_zero, k=2)
        arr[rand[0][0], rand[0][1], rand[0][2]], arr[rand[1][0], rand[1][1], rand[1][2]] = (
            arr[rand[1][0], rand[1][1], rand[1][2]], arr[rand[0][0], rand[0][1], rand[0][2]])
        if new:
            return arr, rand
        else:
            return arr

    def n_pair_swapper(self, n, arr, new=False):
        for i in range(n):
            if new:
                arr, rand = self.pair_swapper(arr, new=new)
                return arr, rand
            else:
                arr = self.pair_swapper(arr, new=new)
                return arr

    def boltzmann_probability(self, delta_e, temperature):
        return np.exp(-delta_e / (self.kb * temperature))

    def logger(self, temp):

        now = datetime.datetime.now()
        now.isoformat()
        system = '-'.join(list(self.initial_config.ele_list))

        folder_path = f'./dump/{system}_{self.initial_config.atoms}_{self.log["ham"]}_{self.config_dict["nn"]}nn'
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
            os.mkdir(f'{folder_path}/plots')

        self.dump_dict = {
            'steps_trajectory': np.array(self.steps),
            'energy_trajectory': np.array(self.energy_trajectory),
            'structure_trajectory': np.array(self.structure_trajectory),
            'Config': self.config_dict,
        }

        plot_conf = {
            "title": f"Energy trajectory for {system} at {temp} K",
            "xlabel": "Steps",
            "ylabel": "Energy (eV/atom)",
            "fontsize": 14,
            "fig_size": (8, 6),
            "file_path": f'./{folder_path}/plots/{temp}_energy.png'
        }
        line_plot(plot_conf=plot_conf, x=np.array(self.steps), y=np.array(self.energy_trajectory))
        #add bar plot

        with open(f'./{folder_path}/{temp}K_{now}.pickle', 'wb') as handle:
            pickle.dump(self.dump_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def mc_single_temp(self, n_trails, temp, lattice, log, nn):
        x_i = lattice.copy()

        if nn == 1:
            e_i = self.energy_finder_new(x_i, self.neighbour_list, flag=log["ham"])
        elif nn == 2:
            e_i = self.energy_finder_new(x_i, self.neighbour_list, flag=log["ham"]) + \
                  (np.sqrt(3) / 2) * self.energy_finder_new(x_i, self.neighbour_list2, flag=log["ham"])

        count = 0
        swaps = 1
        self.steps = []
        self.energy_trajectory = []
        self.structure_trajectory = []
        for i in tqdm(range(n_trails), desc=f"Running at {temp} K with {swaps} swap"):
            x_iplus1, rand = self.n_pair_swapper(n=swaps, arr=x_i.copy(), new=True)
            try:
                assert np.array_equal(x_iplus1, x_i) is not True
            except AssertionError:
                x_iplus1, rand = self.n_pair_swapper(n=1, arr=x_i.copy(), new=True)

            if nn == 1:
                x_i_neighbour = list(set(
                    self.neighbour_dict[tuple(rand[0])] + [tuple(rand[0])] + self.neighbour_dict[tuple(rand[1])] + [
                        tuple(rand[1])]))
                neighbour = [[tuple(i), self.neighbour_dict[i]] for i in x_i_neighbour]
                e_iplus1 = e_i - self.energy_finder_new(x_i, neighbour, flag=log["ham"]) + self.energy_finder_new(
                    x_iplus1.copy(), neighbour, flag=log["ham"])
                e_iplus1 = e_iplus1
            elif nn == 2:
                x_i_neighbour = list(set(
                    self.neighbour_dict[tuple(rand[0])] + [tuple(rand[0])] + self.neighbour_dict[tuple(rand[1])] + [
                        tuple(rand[1])]))
                neighbour = [[tuple(i), self.neighbour_dict[i]] for i in x_i_neighbour]
                x_i_neighbour2 = list(set(
                    self.neighbour_dict2[tuple(rand[0])] + [tuple(rand[0])] + self.neighbour_dict2[tuple(rand[1])] + [
                        tuple(rand[1])]))
                neighbour2 = [[tuple(i), self.neighbour_dict2[i]] for i in x_i_neighbour2]

                e_iplus1 = e_i - self.energy_finder_new(x_i, neighbour, flag=log["ham"]) - (np.sqrt(3) / 2) \
                           * self.energy_finder_new(x_i, neighbour2, flag=log["ham"]) + self.energy_finder_new(
                    x_iplus1.copy(), neighbour, flag=log["ham"]) + (np.sqrt(3) / 2) * self.energy_finder_new(
                    x_iplus1.copy(), neighbour2, flag=log["ham"])
                e_iplus1 = e_iplus1

            if e_iplus1 < e_i or self.boltzmann_probability(e_iplus1 - e_i, temp) >= random.random():
                x_i = x_iplus1
                e_i = e_iplus1

                if count % 10000 == 0:
                    self.steps.append(i)
                    self.energy_trajectory.append(e_iplus1)
                    self.structure_trajectory.append(x_i)

                count += 1

        return x_i

    @property
    def single_temp_protocol(self) -> np.array:
        # system = '-'.join(list(self.initial_config.ele_list))

        print("Warmup Run:")
        x_warm = self.mc_single_temp(n_trails=self.config_dict["n_warm"],
                                     temp=self.config_dict["warm_T"],
                                     lattice=self.initial_config.final_bcclattice,
                                     log=self.log)
        self.steps.append(100)
        self.energy_trajectory.append(100)
        self.structure_trajectory.append(np.zeros_like(x_warm))

        x_final = self.mc_single_temp(n_trails=self.config_dict["n_trails"],
                                      temp=self.config_dict["T"],
                                      lattice=x_warm,
                                      log=self.log)
        self.logger(self.config_dict["T"])

        return x_final

    @property
    def stage_temp_protocol(self) -> np.array:
        temp_ranges = np.linspace(start=3000, stop=300, num=9).astype(int)
        print("Warmup Run:")
        x_warm = self.mc_single_temp(n_trails=self.config_dict["n_warm"],
                                     temp=self.config_dict["warm_T"],
                                     lattice=self.initial_config.final_bcclattice,
                                     log=self.log,
                                     nn=self.config_dict['nn'])
        self.steps.append(100)
        self.energy_trajectory.append(100)
        self.structure_trajectory.append(np.zeros_like(x_warm))

        x_final = x_warm
        for temp in temp_ranges:
            x_final = self.mc_single_temp(n_trails=self.config_dict["n_trails"],
                                          temp=temp,
                                          lattice=x_final,
                                          log=self.log,
                                          nn=self.config_dict['nn'])
            self.logger(temp)

        return x_final


def main(ele_dict: dict,
         ele_assign: dict,
         config_dict: dict,
         lookup: dict,
         rep_unit: int,
         flag: str):
    initial_config = InitialConfig(rep_unit, ele_dict, lookup_dict=lookup, ele_assign=ele_assign)
    print(f"Creating Initial Configuration of {'-'.join(list(ele_dict.keys()))}")
    print(f"{initial_config.atoms} Atoms in total")

    mc = MonteCarlo(initial_config, config_dict)
    if flag == " single_temp":
        _ = mc.single_temp_protocol
    elif flag == "stage_temp":
        _ = mc.stage_temp_protocol


if __name__ == "__main__":
    lookup, ele_assign = enthalpy_model_lookup(source="pravan",
                                               lattice="bcc",
                                               folder_path="/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/data/input_data")
    rep_unit = 10
    ele_dict = {
        'Cr': 0.5,
        'W': 0.5,
    }

    config_dict = {'n_warm': 250000,
                   'warm_T': 4000,
                   'T': 0,
                   'n_trails': 2000000,
                   'log': {'ham': 'bonds'},
                   'nn': 2}

    main(ele_dict=ele_dict,
         config_dict=config_dict,
         lookup=lookup,
         rep_unit=rep_unit,
         flag="stage_temp",
         ele_assign=ele_assign)
