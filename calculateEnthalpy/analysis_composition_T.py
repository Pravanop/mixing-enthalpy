import time
import json
import numpy as np
import pandas as pd
import math
import itertools
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry, Element
from mp_api.client import MPRester
from pymatgen.core import Composition
import matplotlib.pyplot as plt
import seaborn as sns
from extract_intermetallics import extract_intermetallics

"""
json file structure:
comp{
    mol_ratios{
        el1:
        el2:
        ...
    intermetallics{
        ...
    enthalpies{
        ...

"""

def run(el1, el2, el3, add_el):
    df_sqs_2 = pd.read_excel('enthalpy_data_and_predictions/pairwise_mixing_enthalpy.xlsx',
                             sheet_name="our work",
                             index_col=0,
                             nrows=18)
    dfca = pd.DataFrame()
    els = ["Mo", "Nb", "Ta", "V", "W", "Zr", "Ti", "Al", "Hf", "Cr", "C", "Re", "Ru", "Os", "Rh", "Ir", "Si"]
    # els = list(dfca.index)
    comp = []
    for e in els:
        if e in [el1, el2, el3, add_el]:
            comp.append(e)

    # create list of elements in order with respect to els
    # redo into list comprehension

    def pick_two(l):
        out = []
        for i in range(len(l) - 1):
            for j in range(len(l) - i - 1):
                out.append((l[i], l[j + i + 1]))
        return out

    # replace with itertools combinatorics
    dfca["x"] = ""
    dfca["comp"] = ""
    dfca["e1"] = ""
    dfca["e2"] = ""
    dfca["e3"] = ""
    dfca["e4"] = ""
    for i in range(11):
        dfca.loc[i, "x"] = i / 10.0
        dfca.loc[i, "e1"] = comp[0]
        dfca.loc[i, "e2"] = comp[1]
        dfca.loc[i, "e3"] = comp[2]
        dfca.loc[i, "e4"] = comp[3]
        dfca.loc[i, el1] = 1.0 / (3 + i / 10.0)
        dfca.loc[i, el2] = 1.0 / (3 + i / 10.0)
        dfca.loc[i, el3] = 1.0 / (3 + i / 10.0)
        dfca.loc[i, add_el] = (i / 10.0) / (3 + i / 10.0)
        dfca.loc[i, "comp"] = el1 + el2 + el3
        if i != 10:
            dfca.loc[i, "comp"] = dfca.loc[i, "comp"] + add_el + str(i / 10.0)
        else:
            dfca.loc[i, "comp"] = dfca.loc[i, "comp"] + add_el
        dfca.loc[i, "enthalpy"] = 0
        dfca.loc[i, "entropy"] = 0
        # enthalpy
        for pair in pick_two(comp):
            a, b = pair
            if str(df_sqs_2.loc[a, b]) == "nan":  # fix this to be not weird
                dfca.loc[i, "enthalpy"] = dfca.loc[i, "enthalpy"] + df_sqs_2.loc[b, a] * dfca.loc[i, a] * dfca.loc[
                    i, b] * 4.0
            else:
                dfca.loc[i, "enthalpy"] = dfca.loc[i, "enthalpy"] + df_sqs_2.loc[a, b] * dfca.loc[i, a] * dfca.loc[
                    i, b] * 4.0
            # can get rid of if/else statement here, see QuaternaryCombinatorics
        # entropy
        for element in comp:
            molar_fraction = dfca.loc[i, element]
            if molar_fraction != 0:
                dfca.loc[i, "entropy"] = dfca.loc[i, "entropy"] - 8.617333262e-05 * molar_fraction * math.log(
                    molar_fraction)
    dfca.to_excel("4out.xlsx")

    # dfca = pd.read_excel("2out.xlsx")

    def binary(el):
        """
        Return a sorted list of chemical systems
            of the form [...,"Li-Be",...,"Al-Zr",...]
        """
        first_el = el
        second_el = el
        return sorted(["{}-{}".format(*sorted(pair))
                       for pair in itertools.product(first_el, second_el)])

    def ternary(el):
        """
        Return a sorted list of chemical systems
            of the form [...,"Li-Be-S",...,"Al-Zr-Se",...]
        """
        first_el = el
        second_el = el
        third_el = el
        return sorted(["{}-{}-{}".format(*sorted(triple))
                       for triple in itertools.product(first_el, second_el, third_el)])

    # rearrange Materials project database entries into dictionary
    dfim = pd.read_excel("intermetallics_zhaohan.xlsx")
    # dfim.set_index("Unnamed:0", inplace=True)
    im = {}
    for index, row in dfim.iterrows():
        if row["comptype"] not in im:
            im[row["comptype"]] = [(row["comp"], row["Hf"])]
            # Hf -> enthalpy, differentiate DFT enthalpy and mixing_enthalpy
        else:
            im[row["comptype"]].append((row["comp"], row["Hf"]))
    df_3 = pd.read_excel('3z.xlsx')  # see screen_ternary
    df_3.reset_index()
    for index, row in dfca.iterrows():
        e1 = row['e1']
        e2 = row['e2']
        e3 = row['e3']
        e4 = row['e4']
        df_pd = pd.DataFrame()
        # quaternary
        df_pd.at[0, 'comp'] = row['comp']
        df_pd.at[0, 'Hf'] = row['enthalpy']
        df_pd.at[0, 'S'] = row['entropy']
        # ternary
        df_a = df_3[df_3['comp'].str.match(e1 + e2 + e3)]
        df_a.reset_index(drop=True, inplace=True)
        df_b = df_3[df_3['comp'].str.match(e1 + e2 + e4)]
        df_b.reset_index(drop=True, inplace=True)
        df_c = df_3[df_3['comp'].str.match(e1 + e3 + e4)]
        df_c.reset_index(drop=True, inplace=True)
        df_d = df_3[df_3['comp'].str.match(e2 + e3 + e4)]
        df_d.reset_index(drop=True, inplace=True)
        df_pd.at[1, 'comp'] = df_a.loc[0, 'comp']
        df_pd.at[1, 'Hf'] = df_a.loc[0, 'enthalpy']
        df_pd.at[1, 'S'] = df_a.loc[0, 'entropy']
        df_pd.at[2, 'comp'] = df_b.loc[0, 'comp']
        df_pd.at[2, 'Hf'] = df_b.loc[0, 'enthalpy']
        df_pd.at[2, 'S'] = df_b.loc[0, 'entropy']
        df_pd.at[3, 'comp'] = df_c.loc[0, 'comp']
        df_pd.at[3, 'Hf'] = df_c.loc[0, 'enthalpy']
        df_pd.at[3, 'S'] = df_c.loc[0, 'entropy']
        df_pd.at[4, 'comp'] = df_d.loc[0, 'comp']
        df_pd.at[4, 'Hf'] = df_d.loc[0, 'enthalpy']
        df_pd.at[4, 'S'] = df_d.loc[0, 'entropy']
        # rewrite into combinatorics
        # binary
        df_pd.at[5, 'comp'] = e1 + e2
        df_pd.at[5, 'Hf'] = df_sqs_2[e1][e2]
        df_pd.at[5, 'S'] = 5.9730802545007364e-05
        df_pd.at[6, 'comp'] = e1 + e3
        df_pd.at[6, 'Hf'] = df_sqs_2[e1][e3]
        df_pd.at[6, 'S'] = 5.9730802545007364e-05
        df_pd.at[7, 'comp'] = e1 + e4
        df_pd.at[7, 'Hf'] = df_sqs_2[e1][e4]
        df_pd.at[7, 'S'] = 5.9730802545007364e-05
        df_pd.at[8, 'comp'] = e2 + e3
        df_pd.at[8, 'Hf'] = df_sqs_2[e2][e3]
        df_pd.at[8, 'S'] = 5.9730802545007364e-05
        df_pd.at[9, 'comp'] = e2 + e4
        df_pd.at[9, 'Hf'] = df_sqs_2[e2][e4]
        df_pd.at[9, 'S'] = 5.9730802545007364e-05
        df_pd.at[10, 'comp'] = e3 + e4
        df_pd.at[10, 'Hf'] = df_sqs_2[e3][e4]
        df_pd.at[10, 'S'] = 5.9730802545007364e-05
        # rewrite constant into a variable
        # elements
        df_pd.at[11, 'comp'] = e1
        df_pd.at[11, 'Hf'] = 0
        df_pd.at[11, 'S'] = 0
        df_pd.at[12, 'comp'] = e2
        df_pd.at[12, 'Hf'] = 0
        df_pd.at[12, 'S'] = 0
        df_pd.at[13, 'comp'] = e3
        df_pd.at[13, 'Hf'] = 0
        df_pd.at[13, 'S'] = 0
        df_pd.at[14, 'comp'] = e4
        df_pd.at[14, 'Hf'] = 0
        df_pd.at[14, 'S'] = 0
        # query intermetallic from the materials project
        el = {e1, e2, e3, e4}
        df_in = pd.DataFrame()
        cnt = 0
        templist = []
        for pair in binary(el):
            p = pair.replace("-", "")
            if p in im and p not in templist:
                for item in im[p]:
                    df_in.at[cnt, "comp"] = item[0]
                    df_in.at[cnt, "Hf"] = item[1]
                    df_in.at[cnt, "S"] = 0
                    cnt += 1
                templist.append(p)
        for pair in ternary(el):
            p = pair.replace("-", "")
            if p in im and p not in templist:
                for item in im[p]:
                    df_in.at[cnt, "comp"] = item[0]
                    df_in.at[cnt, "Hf"] = item[1]
                    df_in.at[cnt, "S"] = 0
                    cnt += 1
                templist.append(p)
        print(templist)
        # df_pd is solid solution entropic alloys from model
        # df_in is intermetallic, entropy = 0 from materials project
        df_final = pd.concat([df_pd, df_in], ignore_index=True)  # check pd.join?
        # double check comps
        comps = df_final['comp']
        # build phase diagram
        for t in range(13):
            temp = 200 * t
            Ef = df_final['Hf'] - temp * df_final['S']
            mg_comp = [None] * len(comps)
            for i in range(len(comps)):
                mg_comp[i] = Composition(comps[i])
                entries3 = [None] * len(mg_comp)
            for i in range(len(mg_comp)):
                entries3[i] = PDEntry(composition=mg_comp[i], energy=Ef[i] * mg_comp[i].num_atoms)
            entries3.append(PDEntry(composition=Composition("NbV2"), energy=-0.059 * 3))
            # comes from machine learning paper, doi = ??
            print(entries3)
            phase = PhaseDiagram(entries3)
            test = PDEntry(composition=mg_comp[0], energy=Ef[0] * mg_comp[0].num_atoms)
            # test can just be entries3[0]
            dfca.at[index, temp] = float(phase.get_e_above_hull(test)) * 1000
    dfca.to_excel(el1 + el2 + el3 + add_el + ".xlsx")
    # plot using seaform
    data = pd.read_excel(el1 + el2 + el3 + add_el + ".xlsx",
                         usecols="B,N:Z",
                         index_col="x")
    sns.set()
    cmap = sns.cubehelix_palette(start=.5, rot=-.61, light=.98, dark=.35, hue=1, as_cmap=True)
    yticksflip = list(reversed(data.columns))
    data_array = np.flipud(data.transpose())
    sns.heatmap(data_array,
                square=True,
                cmap=cmap,
                xticklabels=data.index,
                yticklabels=yticksflip)
    plt.title("E_hull vs T for " + el1 + el2 + el3 + "+" + add_el + "$_{x}$")
    plt.ylabel("T (K)")
    plt.xlabel(el1 + el2 + el3 + add_el + "$_{x}$")
    plt.ylabel('$E_{hull}$meV/atom)', fontsize=12, labelpad=-30)

    # Add black dashed line where ehull becomes 0
    rows_with_transition = np.where(data_array != 0, 1, 0).sum(axis=0)
    cols_with_transition = np.where(data_array != 0, 1, 0).sum(axis=1)
    # print(rows_with_transition)
    # print(cols_with_transition)
    temp_row = rows_with_transition[0]
    temp_col = cols_with_transition[0]
    for i in range(len(data_array) - 2):
        if temp_row != rows_with_transition[i]:
            height = len(data_array) - temp_row
            plt.plot([temp_col, i], [height, height],
                     color='black',
                     linestyle=(0, (5, 3)),
                     linewidth=1.35)
            plt.plot([i, i], [height, len(data_array) - rows_with_transition[i]],
                     color='black',
                     linestyle=(0, (5, 3)),
                     linewidth=1.35)
            temp_row = rows_with_transition[i]
            temp_col = i
        elif i == len(data_array) - 3:
            height = len(data_array) - temp_row
            plt.plot([temp_col, i + 1], [height, height],
                     color='black',
                     linestyle=(0, (5, 3)),
                     linewidth=1.35)

    plt.show()


run("Nb", "V", "Zr", "Mo")
run("Nb", "V", "Mo", "Zr")
run("Nb", "Zr", "Mo", "V")
run("V", "Zr", "Mo", "Nb")
