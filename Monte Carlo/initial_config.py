import numpy as np


class InitialConfig:
    def __init__(self, rep_unit, ele_dict, lookup_dict):
        self.bcclattice = None
        self.rep_unit = rep_unit
        self.ele_dict = ele_dict
        self.lookup = lookup_dict
        assert sum(list(self.ele_dict.values())) == 1
        self.bcc_lattice()
        self.total_points = self.bcclattice.size
        self.lattice_shape = self.bcclattice.shape
        self.lattice_points = np.where(self.bcclattice == 1)
        self.atoms = np.count_nonzero(self.bcclattice)
        self.ele_list = list(self.ele_dict.keys())
        self.ele_assign = {
            'Cr': 1,
            'W': 2,
            'V': 3,
            'Ta': 4,
            'Ti': 5,
            'Hf': 6,
        }
        self.inv_ele_assign = {v: k for k, v in self.ele_assign.items()}
        self.ele_no_list = [self.ele_assign[i] for i in self.ele_list]
        self.ele_lattice()
        self.indices = np.indices(dimensions=self.lattice_shape).reshape(3, -1).T
        self.non_zero = np.transpose(np.nonzero(self.final_bcclattice))

        self.lookup_reassign()

    def bcc_lattice(self):
        oned_tile1 = [1, 0]
        oned_tile2 = [0, 0]
        oned_tile3 = [0, 1]
        twod_tile1 = [oned_tile1, oned_tile2]
        twod_tile2 = [oned_tile2, oned_tile3]
        threed_tile = [twod_tile1, twod_tile2]
        self.bcclattice = np.tile(threed_tile, (self.rep_unit, self.rep_unit, self.rep_unit))

    def ele_lattice(self):
        elelattice = np.ones(self.total_points).astype(int)

        partition_list = np.round(np.array(list(self.ele_dict.values())) * self.total_points, 0).astype(int)
        chunks = np.split(elelattice, np.cumsum(partition_list[:-1]))
        elelattice_assigned = []
        for idx, chunk in enumerate(chunks):
            elelattice_assigned += list(chunk * self.ele_no_list[idx])
        elelattice_assigned = np.array(elelattice_assigned)
        np.random.shuffle(elelattice_assigned)
        elelattice_assigned = elelattice_assigned.reshape(self.bcclattice.shape)

        self.final_bcclattice = self.bcclattice * elelattice_assigned

    def lookup_reassign(self):
        new_lookup = {}
        for keys, values in self.lookup.items():
            ele_list = keys.split('-')
            revised_key = str(sorted([self.ele_assign[i] for i in ele_list]))
            new_lookup[revised_key] = values

        self.lookup = new_lookup
