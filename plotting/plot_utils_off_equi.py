import pickle
from calculateEnthalpy.helper_functions.grid_code import create_multinary
import matplotlib.pyplot as plt
import matplotlib
from pymatgen.analysis.phase_diagram import PDPlotter

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

			"../data/output_data/old_data/TD_john_from_bokas_offequi_bcc.pickle",
			# "../data/output_data/old_data/TD_bokas_bcc.pickle",
			'rb'
			) as f :
		dump_dict = pickle.load(f)
	phase_diagram = dump_dict[comp]
	print(phase_diagram[1001].stable_entries)
	PDPlotter(phase_diagram[601], show_unstable=True).show()
	# plt.show()



ele_list = ["Cr", "V" , "W",]

compositions = create_multinary(element_list = ele_list , no_comb = [3])[0]
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
