import json
import pickle
from create_alloy_comp import create_multinary
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.core import Composition

font = {
		'family' : 'Helvetica' ,
		'size'   : 16
		}
matplotlib.rc('font' , **font)

def plot_tempDiagram(comp , ax , color) :
	"""

	:param comp:
	"""
	with open(
			"data/output_data/TD_bokas_bcc.pickle" ,
			'rb'
			) as f :
		dump_dict = pickle.load(f)
	data_old = dump_dict[comp]
	data = {}
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
	size = 15
	ax.plot(
			temperatures , energy_above_hull , linestyle = '-' , color = color , alpha = 0.5 ,
			markeredgecolor = "black" , linewidth = 3.5 , label = "_nolegend_" , zorder = 0
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
			s = decomp_latex_name , fontsize = 15
			)
	comp_new = comp.split('-')
	frac = str(round(1 / len(comp_new) , 2))
	print(comp , temperatures[transition])
	join_str = "_" + "{" + frac + "}"
	comp_new_str = "".join(comp_new)
	ax.text(x = temperatures[transition] + 50 , y = 80 , s = comp_new_str , fontsize = 14)
	save = False
	ax.set_xlim([0 , 3000])
	ax.set_ylim([-10 , 250])
	if save :
		# Save or show the plot
		plt.tight_layout()  # Adjust layout to prevent clipping of labels
		plt.savefig('phase_diagram.png' , dpi = 300)  # Save plot as image


ele_list = ["Cr" , "Fe" , "Ta" , "Hf" , "Mn" , "Ti" , "V" , "W"]

compositions = create_multinary(element_list = ele_list , no_comb = [3])[0]
print(compositions)
# compositions = ["Cr-Fe-Hf-Ta-V", "Cr-Fe-Hf-Ta-W", "Cr-Fe-Mn-Ti-V", "Cr-Fe-Ta-Ti-V", "Cr-Fe-Ta-Ti-W"]
# compositions = [ "Cr-Fe-Ta-V-W", "Cr-Fe-Ti-V-W", "Cr-Hf-Ta-V-W", "Cr-Ta-Ti-V-W", "Fe-Hf-Ta-V-W", "Fe-Ta-Ti-V-W",]
# compositions = ["Cr-Fe-Hf-Ta" , "Cr-Fe-Mn-Ti" , "Cr-Fe-Mn-V" ,"Cr-Fe-Ta-Ti","Cr-Fe-Ta-V"]
# compositions =               ["Cr-Fe-V-W" , "Cr-Hf-Ta-W", "Cr-Fe-Ta-W" , "Cr-Fe-Ti-V" , "Cr-Fe-Ti-W" ,]
# compositions =               ["Cr-Hf-Ta-V"   , "Cr-Hf-V-W" , "Cr-Mn-Ti-V" ,"Cr-Ta-Ti-V","Cr-Ta-Ti-W"]
# compositions =               ["Cr-Ta-V-W" , "Cr-Ti-V-W" , "Fe-Hf-Ta-V" ,"Fe-Hf-Ta-W","Fe-Hf-V-W" ]
# compositions =               [ "Fe-Ti-V-W", "Fe-Ta-V-W", "Fe-Mn-Ti-V" , "Fe-Ta-Ti-V" , "Fe-Ta-Ti-W" ,]
# compositions =                ["Hf-Ta-V-W" , "Ta-Ti-V-W"]
fig , axs = plt.subplots(len(compositions) , 1 , sharey = True , sharex = True , figsize = (10 , 8))
plt.subplots_adjust(wspace = 0 , hspace = 0 , left = 0.09 , bottom = 0.095 , top = 0.92 , right = 0.948)
fig.supylabel("Energy above Hull (meV/atom)")
fig.supxlabel("Temperature (K)")
color_scheme = ["#6699CC" , "#004488" , "#EECC66" , "#994455" , "#997700" , "#EE99AA"]
color_scheme += color_scheme
color_scheme += color_scheme
color_scheme += color_scheme
color_scheme += color_scheme
color_scheme += color_scheme
for idx , ax in enumerate(axs) :
	try :
		plot_tempDiagram(compositions[idx] , ax , color_scheme[idx])
	except :
		continue
plt.show()
