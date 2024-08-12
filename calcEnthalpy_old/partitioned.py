import itertools
import os
import random
import time
import pandas as pd
import xlsxwriter
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry, Element
from mp_api.client import MPRester
from pymatgen.core import Composition

'''This is the most I've gotten at combining all the scripts into one.
    The main functions are condensed into these class methods 
    (but you'll likely want to reorganize them further).'''

class AlloySystem:
    '''Creating a new AlloySystem class requires a pandas dataframe.
    This dataframe consists of a 2D grid of binary enthalpies labeled by elements.
    See'''
    def __init__(self,
                 directory: str,
                 datasheet: pd.DataFrame = None,
                 api_key: str = 'NUNc2qkYfekFR1DkxzhKvBCAMVAgOLoF',
                 ):
        self.frames = None
        self.intermetallics_dict = None
        self.mpr = MPRester(api_key=api_key)
        self.directory = directory.rstrip('/') + '/'
        if not os.path.exists(self.directory):
            os.mkdir(self.directory)
        if datasheet is not None:
            self.datasheet = datasheet
            datasheet.to_excel(directory + 'datasheet.xlsx') # Bug: not running
        elif os.path.isfile(directory + 'datasheet.xlsx'):
            self.datasheet = pd.read_excel(directory + 'datasheet.xlsx', index_col=0)
        else:
            raise FileNotFoundError('Binary interaction datasheet not found.')
        self.el = list(self.datasheet.index)

        for els in self.el:
            try:
                Element(els)
            except ValueError:
                raise ValueError('Index values are not all elements.')

    def open(self, sheet: str):
        return pd.read_excel(self.directory + sheet)

    def save(self, df: pd.DataFrame, name):
        df.to_excel(self.directory + name)
        return df

    def binary(self, els):
        first_el = set(els)
        second_el = set(els)
        return sorted({'{}-{}'.format(*sorted(pair))
                       for pair in itertools.product(first_el, second_el)})

    def ternary(self, els):
        first_el = set(els)
        second_el = set(els)
        third_el = set(els)
        return sorted({'{}-{}-{}'.format(*sorted(triple))
                       for triple in itertools.product(first_el, second_el, third_el)})

    def reorder_comp(self, els: list):
        return [e for e in self.el if e in els]

    def DataRetrieval(self, chemsys):
        fields = ['formula_pretty', 'formation_energy_per_atom', 'energy_above_hull']
        docs = self.mpr.summary.search(chemsys=chemsys, fields=fields)
        df = pd.DataFrame(columns=fields)
        for row in range(len(docs)):
            df.at[row, 'formula_pretty'] = docs[row].formula_pretty
            df.at[row, 'formation_energy_per_atom'] = docs[row].formation_energy_per_atom
            df.at[row, 'energy_above_hull'] = docs[row].energy_above_hull
        return df

    def create_intermetallics(self):
        df = pd.DataFrame()
        dfinal = pd.DataFrame()
        temp = []
        for comp in zip(self.binary(self.el), self.ternary(self.el)):
            if comp not in temp:
                df = pd.concat([df, self.DataRetrieval(comp)], ignore_index=False)
                temp.append(comp)
        dfinal['comp'] = df['formula_pretty']
        dfinal['enthalpy'] = df['formation_energy_per_atom']
        dfinal.reset_index(inplace=True, drop=True, allow_duplicates=False)
        for index, row in dfinal.iterrows():
            dfinal.at[index, 'comptype'] = \
                f"{''.join([e for e in sorted(self.el) if Element(e) in Composition(row['comp']).elements])}"
        dfinal.to_excel(self.directory + 'intermetallics.xlsx')
        return

    def intermetallics_to_dict(self):
        path = self.directory + 'intermetallics.xlsx'
        if not os.path.isfile(path):
            print('Retrieving intermetallics from Materials project.')
            self.create_intermetallics()
        dfim = pd.read_excel(path)
        im = {}
        for index, row in dfim.iterrows():
            if row['comptype'] not in im:
                im[row['comptype']] = [(row['comp'], row['enthalpy'])]
            else:
                im[row['comptype']].append((row['comp'], row['enthalpy']))
        return im

    def prep(self, n):
        self.intermetallics_dict = self.intermetallics_to_dict()
        if not os.path.isfile(self.directory + str(n) + "element.xlsx"):
            self.screen(n)
        for nsheet in range(n, 1, -1):
            if not os.path.exists(self.directory + str(nsheet) + 'element.xlsx'):
                print('Creating ' + str(nsheet) + 'element.xlsx')
                self.screen(nsheet)
        self.frames = [self.open(str(n) + 'element.xlsx') for n in range(n, 1, -1)]
        return

    def screen(self, n):
        wb = xlsxwriter.Workbook(self.directory + str(n) + 'element.xlsx',
                                 {'strings_to_numbers': True})
        length = len(self.datasheet.index)
        s1 = wb.add_worksheet()
        s1.write(0, 0, 'comp')
        for i in range(n + 1):
            s1.write(0, i + 1, 'e' + str(i + 1))
        s1.write(0, n + 1, 'enthalpy')
        s1.write(0, n + 2, 'entropy')
        index = 1
        combinations = itertools.combinations(range(length), n)
        for c in combinations:
            comp = ''
            var = []
            for i in range(n):
                s1.write(index, i + 1, self.datasheet.columns.values[c[i]])
                var.append(self.datasheet.columns.values[c[i]])
                comp += self.datasheet.columns.values[c[i]]
            s1.write(index, 0, comp)
            enthalpy_sum = 0
            for i in range(n):
                for j in range(i + 1, n):
                    add = self.datasheet.loc[var[j], var[i]]
                    enthalpy_sum += float(add)
            s1.write(index, n + 1, enthalpy_sum / (n ** 2) * 4)
            GAS_CONSTANT_R = 8.617333262e-05
            s1.write(index, n + 2, GAS_CONSTANT_R * math.log(n))
            index += 1
        wb.close()
        return str(n) + 'element.xlsx'

    def create_df_pd(self, elements: list, intermetallics=True):
        print(elements)
        elements = self.reorder_comp(elements)  # Check: may not be necessary
        alloy_size = len(elements)
        df_pd = pd.DataFrame()
        pd_index = 0
        for i in range(alloy_size, 0, -1):
            if i != 1:
                df = self.frames[alloy_size - i]
                element_combinations = itertools.combinations(elements, i)
                for comp in list(element_combinations):
                    dfn = df[df['comp'].str.fullmatch(''.join(comp))] # str.match > str.fullmatch (edge-case)
                    dfn.reset_index(drop=True, inplace=True)
                    df_pd.at[pd_index, 'comp'] = dfn.loc[0, 'comp']
                    df_pd.at[pd_index, 'enthalpy'] = dfn.loc[0, 'enthalpy']
                    df_pd.at[pd_index, 'entropy'] = dfn.loc[0, 'entropy']
                    pd_index += 1
            else:
                for el_num in range(alloy_size):
                    df_pd.at[pd_index, 'comp'] = elements[el_num]
                    df_pd.at[pd_index, 'enthalpy'] = 0
                    df_pd.at[pd_index, 'entropy'] = 0
                    pd_index += 1
        if intermetallics:
            for comp in itertools.chain(self.binary(elements), self.ternary(elements)):
                comp = comp.replace('-', '')
                if comp in self.intermetallics_dict:
                    for item in self.intermetallics_dict[comp]:
                        df_pd.at[pd_index, 'comp'] = item[0]
                        df_pd.at[pd_index, 'enthalpy'] = item[1]
                        df_pd.at[pd_index, 'entropy'] = 0
                        pd_index += 1
        return df_pd

    def phase_diagram(self, df_final, temperature):
        comps = df_final['comp']
        formation_energy = df_final['enthalpy'] - temperature * df_final['entropy']
        mg_comp = [None] * len(comps)
        for i in range(len(comps)):
            mg_comp[i] = Composition(comps[i])
            entries3 = [None] * len(mg_comp)
        for i in range(len(mg_comp)):
            entries3[i] = PDEntry(composition=mg_comp[i], energy=formation_energy[i] * mg_comp[i].num_atoms)
        phase = PhaseDiagram(entries3)
        test = PDEntry(composition=mg_comp[0], energy=formation_energy[0] * mg_comp[0].num_atoms)
        decomp, ehull = phase.get_decomp_and_e_above_hull(test)
        return float(ehull), decomp

    def run_ehull_combinations(self, n, temperature=0):
        self.prep(n)
        df = self.open(str(n) + "element.xlsx")
        df['phase'] = ""
        df['ehull'] = ""
        for index, row in df.iterrows():
            current_elements_list = [row['e' + str(i + 1)] for i in range(n)]
            df_pd = self.create_df_pd(current_elements_list)
            energy_above_hull, phase = self.phase_diagram(df_pd, temperature)
            df.at[index, 'phase'] = phase
            df.at[index, 'ehull'] = energy_above_hull
        self.save(df, str(n) + '-e' + str(temperature) + '.xlsx')
        return

    def find_stable_temp(self, ):
        temperature = []
        eahull = []
        decomps = []
        for temp in range(0, 2050, 50):
            temperature.append(temp)
            x, phase = phasecheck(temp)
            decomps.append((phase, temp))
            eahull.append(x)
        fig, ax = plt.subplots()
        ax.plot(temperature, eahull)
        ax.fill_between(temperature, eahull, 0, alpha=0.3)
        plt.xlabel('Temperature (K)')
        plt.ylabel('Energy above hull')
        tempphase = ''
        c = 0
        for d in decomps:
            if decompsplit(d[0]) != tempphase:
                ax.text(d[1], 0.001 * c, decompsplit(d[0]))
                tempphase = decompsplit(d[0])
                c += 1
        title = 'Stability vs. Temperature for ' + el1 + el2 + el3 + el4 + el5
        plt.title(title)
        plt.show()
        return


def decompsplit(string):
    begin = []
    end = []
    comps = []
    string = string.replace('PDEntry : ', '|')
    string = string.replace(' with energy ', '?')
    for i in range(len(string)):
        if string[i] == '|':
            begin.append(i)
        if string[i] == '?':
            end.append(i)
    for i in range(len(begin)):
        comps.append(Composition(string[begin[i] + 1: end[i]]).reduced_formula)
    return comps
