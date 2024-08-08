import itertools
import time
import pandas as pd
import math
from pymatgen.core import Composition
from pymatgen.ext.matproj import MPRester
import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry, Element
from matminer.data_retrieval.retrieve_MP import MPDataRetrieval
from pymatgen.entries.compatibility import MaterialsProjectCompatibility

# screen_binary.py: Creates base document for pairwise enthalpies
def r(structure):
    df_sqs_2 = pd.read_excel('bokas.xlsx', sheet_name=structure, index_col=0)
    els = ['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn', 'Mo',
           'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr', 'V', 'Mg', 'Re',
           'Os', 'Rh', 'Ir', 'Pd', 'Pt', 'Ag', 'Au', 'Zn', 'Cd']
    arr = pd.DataFrame(columns=["comp", "e1", "e2", "enthalpy", "entropy"])
    r = 0
    for a in els:
        for b in els:
            if a != b and els.index(a) < els.index(b):
                arr.at[r, "comp"] = a + b
                arr.at[r, "e1"] = a
                arr.at[r, "e2"] = b
                arr.at[r, "enthalpy"] = df_sqs_2.at[b, a]
                arr.at[r, "entropy"] = 0.0000597
                r += 1
    arr.to_excel("2element" + structure + ".xlsx", sheet_name=structure)


# r("FCC")
# r("BCC")
# sheets need to be combined into 2element.xlsx
# update2.py: Creates phase diagrams and determines e_hull + decompositions
def run(structure):
    df_sqs_2 = pd.read_excel('bokas.xlsx', sheet_name=structure)
    df_sqs_2.set_index('Unnamed: 0', inplace=True)

    def binary():
        """
        Return a sorted list of chemical systems
            of the form [...,"Li-Be-S",...,"Al-Zr-Se",...]
        """
        first_el = el
        second_el = el
        return sorted(["{}-{}".format(*sorted(pair))
                       for pair in itertools.product(first_el, second_el)])

    df_2 = pd.read_excel('2element.xlsx', sheet_name=structure)
    for index, row in df_2.iterrows():
        e1 = row['e1']
        e2 = row['e2']
        df_pd = pd.DataFrame()
        # binary
        df_pd.at[0, 'comp'] = row['comp']
        df_pd.at[0, 'Hf'] = row['enthalpy']
        df_pd.at[0, 'S'] = row['entropy']
        # elements
        df_pd.at[1, 'comp'] = e1
        df_pd.at[1, 'Hf'] = 0
        df_pd.at[1, 'S'] = 0
        df_pd.at[2, 'comp'] = e2
        df_pd.at[2, 'Hf'] = 0
        df_pd.at[2, 'S'] = 0
        # query binary intermetallic from the materials project
        el = {e1, e2}
        dbin = pd.read_excel("database_binary.xlsx")
        dbin.set_index("index", inplace=True)
        df_in = pd.DataFrame()
        cnt = 0
        for pair in binary():
            a, b = pair.split("-")
            if a != b:
                for i, r in dbin.iterrows():
                    if r["comp"].find(a) != -1 and r["comp"].find(b) != -1:
                        df_in.at[cnt, "comp"] = r["comp"]
                        df_in.at[cnt, "Hf"] = r["Hf"]
                        df_in.at[cnt, "S"] = r["S"]
                        cnt += 1
        df_final = pd.concat([df_pd, df_in], ignore_index=True)
        # build phase diagram
        comps = df_final['comp']
        Ef = df_final['Hf'] - 1000 * df_final['S']
        mg_comp = [None] * len(comps)
        for i in range(len(comps)):
            mg_comp[i] = Composition(comps[i])
            entries3 = [None] * len(mg_comp)
        for i in range(len(mg_comp)):
            entries3[i] = PDEntry(composition=mg_comp[i], energy=Ef[i])
        phase = PhaseDiagram(entries3)
        # get decomposition and e_above_hull
        test = PDEntry(composition=mg_comp[0], energy=Ef[0])
        out = str(phase.get_decomp_and_e_above_hull(test))
        x = out.split("}, ")[1].rstrip(")")
        df_2.at[index, "e_hull"] = x
        y = out.split("}, ")[0].lstrip("(") + "}"
        print(y)
        df_2.at[index, "phase"] = y
    df_2.to_excel("2-" + structure.lower() + ".xlsx")
    return -1


run("FCC")
run("BCC")
