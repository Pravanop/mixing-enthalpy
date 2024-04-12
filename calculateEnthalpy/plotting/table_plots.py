import pickle

import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pymatgen.core import Composition
from prepareVASPRuns.file_utils import load_json_to_dict

font = {
		'family' : 'Helvetica' ,
		'size'   : 16
		}
matplotlib.rc('font' , **font)
file_path = "../data/output_data/old_data/binary_bcc_Cr_W_equiOut.csv"
df = pd.read_csv(file_path)
# df['mix_enthalpy (meV/atom)'] = round(df['mix_enthalpy (meV/atom)']*1000,0)
single_energy = load_json_to_dict("../data/input_data/single_energy.json")

comp_list = np.array(df["Comp"])
avg_tm = []
misc_t = []
decomp_list = []
for comp in comp_list:
	ele_list = comp.split('-')
	temp = 0
	for j in ele_list:
		temp += single_energy[j]['Tm']
	
	temp /= len(ele_list)
	avg_tm.append(round(temp, 0))

	# with open(
	# 		"data/output_data/TD_bokas_bcc.pickle" ,
	# 		'rb'
	# 		) as f :
	# 	dump_dict = pickle.load(f)
	# data_old = dump_dict[comp]
	# data = {}
	# for key , value in data_old.items() :
	# 	data[key] = value[0] * 1000
	# decomp = {}
	# enthalpy = {}
	# for key , value in data_old.items() :
	# 	decomp[key] = value[1]
	# 	enthalpy[key] = value[2]
	#
	# decomp_temp = list(decomp.values())
	# temperatures = np.array(list(data.keys())).astype(int)
	# energy_above_hull = list(data.values())
	# try :
	# 	transition = energy_above_hull.index(0)
	# except :
	# 	transition = -1
	# decomp_transition = decomp_temp[transition - 1]
	# decomp_latex_name = ""
	# for key in list(decomp_transition.keys()) :
	# 	decomp_latex_name += Composition(
	# 			''.join(
	# 					str(key).split(": ")[1].split(" w")[0].split(
	# 							" "
	# 							""
	# 							)
	# 					)
	# 			).reduced_formula + " + "
	#
	# decomp_latex_name = decomp_latex_name[:-3]
	# misc_t.append(temperatures[transition] - 1)
	# decomp_list.append(decomp_latex_name)

df['Tm'] = avg_tm
# df['Misc_T'] = misc_t
# df['Decomposition'] = decomp_list
df_proc = df.loc[df['Tmisc'] <= 20*df['Tm']]
df_proc = df_proc.sort_values(['Tm'])
# plt.grid(axis = 'y', zorder = 0)
fig, ax = plt.subplots()
df_proc.plot.bar(x='Comp', ax=ax, color = ['#004488','#994455','#6699CC' ], edgecolor = 'black')
print(df_proc.head())
# plt.bar(x = df_proc['comp'], height = df_proc['Tm'], color = '#004488', zorder = 3, width = 0.2, edgecolor = 'black')
# plt.bar(x = df_proc['comp'], height = df_proc['Misc_T'], color = '#994455', zorder = 3, width = 0.2, edgecolor =
# 'black')
# plt.bar(x = df_proc['Comp'], height = df_proc['Expt'], color = '#6699CC', zorder = 3, width = 0.2)
plt.xticks(rotation = 45)
plt.ylabel("Temperature (K)")

plt.xlabel("Composition")
plt.legend([ "$T_{misc}$", "Expt.", "$Avg. T_{m}$",])
plt.show()


# df.drop(columns = "config_entropy", inplace = True)
# df.to_csv(file_path[:-4] + "_processed.csv")