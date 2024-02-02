import itertools
import random
import time
import pandas as pd
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry, Element
from matminer.data_retrieval.retrieve_MP import MPDataRetrieval
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.core import Composition
from mp_api.client import MPRester
import plotly.graph_objects as go
import re
#
# df = pd.DataFrame()
# for x in range(4):
#     number = x + 2  # 2, 3, 4, 5
#     for i in range(11):
#         y = i / 10.0
#         element = 1 / (number + y)
#         dependent = y / (number + y)
#         S = number * element * math.log(element) + dependent * math.log(element)
#         df.loc[x, y] = S
# df.to_excel("1o.xlsx")

directory = ""

dfim = pd.read_excel(directory + "intermetallics.xlsx")
dfim.set_index("index", inplace=True)
im = {}
for index, row in dfim.iterrows():
    if row["comptype"] not in im:
        im[row["comptype"]] = [(row["comp"], row["Hf"])]
    else:
        im[row["comptype"]].append((row["comp"], row["Hf"]))

def binary():
    """
    Return a sorted list of chemical systems
        of the form [...,"Li-Be-S",...,"Al-Zr-Se",...]
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

def decompsplit(string):
    begin = []
    end = []
    comps = []
    string = string.replace("PDEntry : ", "|")
    string = string.replace(" with energy ", "?")
    for i in range(len(string)):
        if string[i] == "|":
            begin.append(i)
        if string[i] == "?":
            end.append(i)
    for i in range(len(begin)):
        comps.append(Composition(string[begin[i] + 1: end[i]]).reduced_formula)
    return comps

structure = "FCC"
df_sqs_2 = pd.read_excel('bokas.xlsx', sheet_name=structure)
df_sqs_2.set_index('Unnamed: 0', inplace=True)

df_3 = pd.read_excel(directory + '3element.xlsx', sheet_name='Ternary ' + structure)
df_4 = pd.read_excel(directory + '4element.xlsx', sheet_name='Quaternary ' + structure)
df_4["phase"] = ""
el1 = "Cu"
el2 = "Ag"
el3 = "Au"
el4 = "Zn"
for index, row in df_4.iterrows():
    e1 = row['e1']
    e2 = row['e2']
    e3 = row['e3']
    e4 = row['e4']
    if str(e1 + e2 + e3 + e4) == str(el1 + el2 + el3 + el4):
        df_pd = pd.DataFrame()
        # quarternary
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
        # query binary intermetallic from the materials project
        el = {e1, e2, e3, e4}
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
        temp = 0

        def phasecheck(t):
            Ef = df_final['Hf'] - t * df_final['S']
            mg_comp = [None] * len(comps)
            for i in range(len(comps)):
                mg_comp[i] = Composition(comps[i])
                entries3 = [None] * len(mg_comp)
            for i in range(len(mg_comp)):
                entries3[i] = PDEntry(composition=mg_comp[i], energy=Ef[i] * mg_comp[i].num_atoms)
            phase = PhaseDiagram(entries3)
            df_4.loc[index, 'phase'] = phase
            # get decomposition and e_above_hull
            test = PDEntry(composition=mg_comp[0], energy=Ef[0] * mg_comp[i].num_atoms)
            entry = str(phase.get_decomp_and_e_above_hull(test))
            ehull = entry.split("}, ")[1].rstrip(")")
            decomp = entry.split("}, ")[0].lstrip("(") + "}"
            df_4.loc[index, "e_hull"] = ehull
            return float(ehull), decomp

        temperature = []
        eahull = []
        decomps = []
        for x in range(50):
            temperature.append(temp)
            x, phase = phasecheck(temp)
            decomps.append((phase, temp))
            eahull.append(x)
            temp += 50
        # fig, ax = plt.subplots()
        # ax.plot(temperature, eahull)
        # ax.fill_between(temperature, eahull, 0, alpha=0.3)
        # plt.xlabel('Temperature (K)')
        # plt.ylabel('Energy above hull')

        # title = "Stability vs. Temperature for " + el1 + el2 + el3
        # plt.title(title)
        # plt.show()

        # Create a line plot
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=temperature,
            y=eahull,
            mode='lines',
            line=dict(
                color='blue',
            ),
            name='Energy Above Hull',
        ))
        tempphase = ""
        for d in decomps:
            if str(decompsplit(d[0])) != tempphase:
                tempphase = str(decompsplit(d[0]))
                label_text = tempphase.replace("[", "").replace("]", "").replace("'", "")
                phase_change_temp = d[1]
                phase_change_energy = eahull[decomps.index(d)]
                fig.add_trace(go.Scatter(
                    x=[phase_change_temp],
                    y=[phase_change_energy],
                    mode='markers',
                    marker=dict(
                        size=10,
                        color='red',
                        symbol='x',
                    ),
                    textposition='top center',
                    name="Phase Change",
                    showlegend=False
                ))
                fig.add_annotation(
                    x=phase_change_temp,
                    y=phase_change_energy,
                    text=label_text,
                    showarrow=True,
                    arrowhead=1,
                    ax=20,
                    ay=-40,
                )

        # Set the layout and labels
        comp = ""
        for e in el:
            comp += e
            comp += " "
        fig.update_layout(
            title='Temperature vs Energy Above Hull for ' + comp,
            xaxis_title='Temperature (K)',
            yaxis_title='Energy Above Hull',
            legend=dict(
                orientation='h',
                yanchor='top',
                y=1.1,
                xanchor='center',
                x=0.5,
            ),
        )

        # Show the graph
        fig.show()

