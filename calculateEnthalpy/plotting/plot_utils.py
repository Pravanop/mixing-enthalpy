import json
import pickle
from calculateEnthalpy.create_alloy_comp import create_multinary
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.core import Composition

font = {
		'family' : 'Helvetica' ,
		'size'   : 20
		}
matplotlib.rc('font' , **font)

def plot_tempDiagram(comp , ax , color, save) :
	"""

	:param comp:
	"""
	with open(
			"../data/output_data/old_data/TD_bokas_bcc.pickle",
			'rb'
			) as f :
		dump_dict = pickle.load(f)
	data_old = dump_dict[comp]
	data = {}
	print(data_old)
	for key , value in data_old.items() :
		data[key] = value[0] * 1000
	decomp = {}
	enthalpy = {}
	for key , value in data_old.items() :
		decomp[key] = value[1]
		enthalpy[key] = value[2]
	
	decomp_temp = list(decomp.values())
	temperatures = np.array(list(data.keys())).astype(int)
	energy_above_hull = list(data.values())
	try :
		transition = energy_above_hull.index(0)
	except :
		transition = -1
	decomp_transition = decomp_temp[transition - 1]
	decomp_latex_name = ""
	for key in list(decomp_transition.keys()) :
		decomp_latex_name += Composition(
				''.join(
						str(key).split(": ")[1].split(" w")[0].split(
								" "
								""
								)
						)
				).to_latex_string() + " + \n"
	
	decomp_latex_name = decomp_latex_name[:-3]
	size = 50
	ax.plot(
			temperatures , energy_above_hull , linestyle = '-' , color = color , alpha = 0.5 ,
			markeredgecolor = "black" , linewidth = 5.5 , label = "_nolegend_" , zorder = 0
			)
	ax.scatter(
			temperatures , energy_above_hull , marker = 'o' , linestyle = '-' , color = color , s = size * 9 ,
			edgecolor = 'black'
			)
	ax.axhline(y = 0 , linewidth = 1.5 , linestyle = '--' , color = 'black' , zorder = 0)
	ax.axhline(y = 50 , linewidth = 2 , linestyle = '-' , color = '#882255' , zorder = 0)
	ax.axvline(x = temperatures[transition] , linewidth = 2.2 , linestyle = '--' , color = 'black' , zorder = 0)
	ax.axvspan(0 , temperatures[transition] , alpha = 0.1 , color = color)
	ax.text(
			x = temperatures[transition] - 50 , y = 50 , verticalalignment = "bottom" , horizontalalignment = "right" ,
			s = decomp_latex_name
			)
	comp_new = comp.split('-')
	frac = str(round(1 / len(comp_new) , 2))
	print(comp , temperatures[transition])
	join_str = "_" + "{" + frac + "}"
	comp_new_str = "".join(comp_new)
	ax.text(x = temperatures[transition] + 50 , y = 80 , s = comp_new_str)
	ax.set_xlim([0 , 2000])
	ax.set_ylim([-10 , 150])
	if save :
		# Save or show the plot
		# plt.tight_layout()  # Adjust layout to prevent clipping of labels
		plt.savefig('/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/plotting/phase_diagram.png' , dpi = 300)  # Save plot as image


ele_list = ["Cr" , "Ta" , "V" , "W", "Ti"]

compositions = create_multinary(element_list = ele_list , no_comb = [5])[0]

fig , axs = plt.subplots(len(compositions) , 1 , sharey = True , sharex = True , figsize = (10 , 8))
if len(compositions) == 1:
	axs = [axs]
plt.subplots_adjust(wspace = 0 , hspace = 0 , left = 0.12 , bottom = 0.095 , top = 0.92 , right = 0.948)
fig.supylabel("Energy above Hull (meV/atom)")
fig.supxlabel("Temperature (K)")
color_scheme = ["#6699CC" , "#004488" , "#EECC66" , "#994455" , "#997700" , "#EE99AA"]
for idx , ax in enumerate(axs) :
	# try :
	plot_tempDiagram(compositions[idx] , ax , color_scheme[idx], save = True)
	# except :
	# 	continue
# plt.show()
