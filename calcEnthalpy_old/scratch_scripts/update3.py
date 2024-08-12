import itertools
import pandas as pd
import math
import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry, Element
from mp_api.client import MPRester
from pymatgen.core import Composition


def run(structure):
    df_sqs_2 = pd.read_excel('bokas.xlsx', sheet_name=structure)
    df_sqs_2.set_index('Unnamed: 0', inplace=True)
    # df_sqs_2 = pd.read_excel('enthalpy_data_and_predictions/pairwise_mixing_enthalpy.xlsx',
    #                          sheet_name="our work",
    #                          index_col=0,
    #                          nrows=18)

    def binary():
        """
        Return a sorted list of chemical systems
            of the form [...,"Li-Be",...,"Al-Zr",...]
        """
        first_el = el
        second_el = el
        return sorted(["{}-{}".format(*sorted(pair))
                       for pair in itertools.product(first_el, second_el)])

    def ternary():
        """
        Return a sorted list of chemical systems
            of the form [...,"Li-Be-S",...,"Al-Zr-Se",...]
        """
        first_el = el
        second_el = el
        third_el = el
        return sorted(["{}-{}-{}".format(*sorted(triple))
                       for triple in itertools.product(first_el, second_el, third_el)])

    df_3 = pd.read_excel('3element.xlsx', sheet_name='Ternary ' + structure)
    # df_3 = pd.read_excel('3z.xlsx')
    df_3["phase"] = ""
    # rearrange Materials project database entries into dictionary
    dfim = pd.read_excel("intermetallics.xlsx")
    # dfim = pd.read_excel("intermetallics_zhaohan.xlsx")
    dfim.set_index("index", inplace=True)
    im = {}
    for index, row in dfim.iterrows():
        if row["comptype"] not in im:
            im[row["comptype"]] = [(row["comp"], row["Hf"])]
        else:
            im[row["comptype"]].append((row["comp"], row["Hf"]))
    for index, row in df_3.iterrows():
        e1 = row['e1']
        e2 = row['e2']
        e3 = row['e3']
        df_pd = pd.DataFrame()
        # ternary
        df_pd.at[0, 'comp'] = row['comp']
        df_pd.at[0, 'Hf'] = row['enthalpy']
        df_pd.at[0, 'S'] = row['entropy']
        # binary
        df_pd.at[1, 'comp'] = e1 + e2
        df_pd.at[1, 'Hf'] = df_sqs_2[e1][e2]
        df_pd.at[1, 'S'] = 5.9730802545007364e-05
        df_pd.at[2, 'comp'] = e1 + e3
        df_pd.at[2, 'Hf'] = df_sqs_2[e1][e3]
        df_pd.at[2, 'S'] = 5.9730802545007364e-05
        df_pd.at[3, 'comp'] = e2 + e3
        df_pd.at[3, 'Hf'] = df_sqs_2[e2][e3]
        df_pd.at[3, 'S'] = 5.9730802545007364e-05
        # elements
        df_pd.at[4, 'comp'] = e1
        df_pd.at[4, 'Hf'] = 0
        df_pd.at[4, 'S'] = 0
        df_pd.at[5, 'comp'] = e2
        df_pd.at[5, 'Hf'] = 0
        df_pd.at[5, 'S'] = 0
        df_pd.at[6, 'comp'] = e3
        df_pd.at[6, 'Hf'] = 0
        df_pd.at[6, 'S'] = 0
        # query binary intermetallic from the materials project
        el = {e1, e2, e3}
        df_in = pd.DataFrame()
        cnt = 0
        templist = []
        for pair in binary():
            p = pair.replace("-", "")
            if p not in templist:
                if p in im and p not in templist:
                    for item in im[p]:
                        df_in.at[cnt, "comp"] = item[0]
                        df_in.at[cnt, "Hf"] = item[1]
                        df_in.at[cnt, "S"] = 0
                        cnt += 1
                    templist.append(p)
        for ter in ternary():
            p = ter.replace("-", "")
            if p in im and p not in templist:
                for item in im[p]:
                    df_in.at[cnt, "comp"] = item[0]
                    df_in.at[cnt, "Hf"] = item[1]
                    df_in.at[cnt, "S"] = 0
                    cnt += 1
                templist.append(p)
        df_final = pd.concat([df_pd, df_in], ignore_index=True)
        # build phase diagram
        comps = df_final['comp']
        Ef = df_final['Hf'] - 1000 * df_final['S']
        mg_comp = [None] * len(comps)
        for i in range(len(comps)):
            mg_comp[i] = Composition(comps[i])
            entries3 = [None] * len(mg_comp)
        for i in range(len(mg_comp)):
            entries3[i] = PDEntry(composition=mg_comp[i], energy=Ef[i]*mg_comp[i].num_atoms)
        phase = PhaseDiagram(entries3)
        df_3.loc[index, 'phase'] = phase
        # get decomposition and e_above_hull
        test = PDEntry(composition=mg_comp[0], energy=Ef[0]*mg_comp[0].num_atoms)
        out = str(phase.get_decomp_and_e_above_hull(test))
        x = out.split("}, ")[1].rstrip(")")
        df_3.at[index, "e_hull"] = x
        y = out.split("}, ")[0].lstrip("(") + "}"
        df_3.at[index, "phase"] = y
    df_3.to_excel("3newbcc.xlsx")
    return 1


run("BCC")
