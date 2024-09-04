import numpy as np
import streamlit as st
from pymatgen.core import Element

from UI.search import search_composition
from pymatgen.analysis.phase_diagram import PDPlotter
import pandas as pd

from calculateEnthalpy.UMAPS.use_umaps import use_umaps_Tmisc
from calculateEnthalpy.heatmap_plots import add_ele
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram
from calculateEnthalpy.reactionPathways.new_rP import new_rP


def convert_df(df):
	# IMPORTANT: Cache the conversion to prevent computation on every rerun
	return df.to_csv().encode('utf-8')


def clicked(button):
	st.session_state.clicked[button] = True

def on_change():
	st.session_state.clicked[1] = False

if 'clicked' not in st.session_state:
	st.session_state.clicked = dict(zip(range(-10, 10), [False] * 20))

element_list_path = (
	"/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/data/input_data/bokasCorrected/element_list_BCC_bokasCorrected.txt")

with open(element_list_path) as f:
	element_total_list = f.read()
element_total_list = element_total_list.split(',')

st.set_page_config(page_title="FAR-HEAA",
				   initial_sidebar_state="expanded",
				   layout="wide")

st.title('Fast and Robust High Entropy Alloy Analysis')
len_check = True
rep_check = True
inv_check = True
buff, col, buff2 = st.columns([1, 3, 1])
user_inp = col.text_input(label='input',
						  placeholder="For example: Al-Cu-Ni",
						  label_visibility='hidden')
ele_list_user_inp = user_inp.split('-')

print(ele_list_user_inp)

if ele_list_user_inp != ['']:
	len_check = set(ele_list_user_inp) <= set(element_total_list)

if not len_check:
	st.write(
		f":red["
		f":blue[{', '.join(list(set(ele_list_user_inp).difference(set(element_total_list))))}] not in "
		f"database! We are working on increasing our database!]"
	)
	inv_check = False

rep_check = len(ele_list_user_inp) == len(set(ele_list_user_inp))
if rep_check is False:
	st.write(":red[Repeated Element!]")

one_check = True
if len(ele_list_user_inp) == 1:
	st.write(":red[Please provide more Elements!]")
	one_check = False

invalid = True
col1, col2 = st.columns([1, 1])
if user_inp and rep_check and len_check and inv_check and one_check:

	with col1:
		# find_comp = st.button('Convex Hulls for this Composition', on_click=clicked, args=[2])
		find_heatmap = st.toggle('Pairwise Mixing Enthalpy', args=[3])

		find_comp = st.toggle('Convex Hulls', args=[2])
		find_misc_T = st.toggle('Miscibility Temperature', args=[4])

		addition_ele = st.toggle('Add an element', args=[6])
		find_reaction_pathway = st.toggle('Deposition Pathways', args=[5])
		UMAP_viz = st.toggle('Visualize UMAP', args=[7])
		decomposition_products = st.toggle('Decomp Products', args=[8])

	with col2:

		if find_comp or find_misc_T or find_heatmap or find_reaction_pathway or addition_ele or UMAP_viz or decomposition_products:
			sys_opts = ["Biasing", "Include_IM", "Include_only_equimolar"]

			col_checkbox = st.columns([1, 1, 1.5])
			alloy_sys = {symbol: col_checkbox[i].checkbox(symbol) for i, symbol in enumerate(sys_opts)}
			find = list(alloy_sys.values())

			if UMAP_viz:
				find[2] = True

			if find[0]:
				binary_file_path = "calculateEnthalpy/new_phase_diagram/bokas_omegas_processed.json"
			else:
				binary_file_path = "data/output_data/bokasCorrected_bcc_1/all_lattices_binaries.json"

			end_member_path = "calculateEnthalpy/new_phase_diagram/bokas_end_members_dict.json"

			pD = phaseDiagram(
				processed_binary_file_path=binary_file_path,
				end_member_file_path=end_member_path,
				grid_size=10,
				im_flag=find[1],
				correction=find[0],
				equi_flag=find[2])

		if find_comp:
			T = st.slider(label="Temperature (K)", min_value=0, max_value=3000)
			conv_hull = pD.make_convex_hull(composition=ele_list_user_inp, temperature=T)
			st.write(PDPlotter(conv_hull).get_plot(process_attributes=True))
			series1  = []
			series2 = []
			for i in conv_hull.stable_entries:
				series1.append(i.name)
				series2.append(i.energy_per_atom)
			df = pd.DataFrame([series1, series2]).T
			df.columns = ['Phase', 'Energy/per_atom']
			st.write(df)

		if find_misc_T:

			mol_user_inp = st.text_input(label='Enter Mole Fractions',
										 placeholder="For example: 0.2-0.3-0.5",
										 label_visibility='visible')
			if mol_user_inp != '':
				mol_user_inp_list = np.array(mol_user_inp.split('-')).astype(float)
			else:
				mol_user_inp_list = []

			print(mol_user_inp_list)
			if len(mol_user_inp_list) < len(ele_list_user_inp):
				st.write(":red[Invalid Input]")
			elif len(mol_user_inp_list) <= 1:
				st.write(":red[No Input]")
			elif sum(mol_user_inp_list) != 1:
				st.write(":red[mole fractions must always sum upto 1")
			else:
				misc_T = pD.find_misc_temperature(composition=ele_list_user_inp,
												  mol_ratio=mol_user_inp_list,
												  lattice='min')
				if isinstance(misc_T, float):
					st.write(f"Miscible Temperature: {misc_T} K")
					if misc_T > 200:
						st.write(f"Decomposition products at 200K lesser than {misc_T} K")
						st.write(pD.find_decomp_products(composition=ele_list_user_inp,
														 mol_ratio=mol_user_inp_list,
														 temperature=misc_T-200,
														 lattice='min')[0])

				else:
					st.write(misc_T + "K")

		if addition_ele:
			add_el_user_inp = st.text_input(label='Enter alloying ele',
											placeholder="For example: Al",
											label_visibility='visible')
			if add_el_user_inp != '':
				ax = add_ele(composition=ele_list_user_inp,
						add_el=add_el_user_inp,
						pD=pD)

				st.write(ax.get_figure())

		if find_heatmap:
			st.write(pD.heatmap(ele_list_user_inp).get_figure())

		if UMAP_viz:
			if len(ele_list_user_inp) == 3:
				st.pyplot(use_umaps_Tmisc(composition=ele_list_user_inp,
										 pD=pD,
										 n=len(ele_list_user_inp))[1])
			else:
				st.write(':red[Currently only works for ternaries. Come back later for higher order visualizations.]')

		if find_reaction_pathway:
			ax, fig = new_rP(composition=ele_list_user_inp,pD=pD)
			st.pyplot(fig)

		if decomposition_products:
			mol_user_inp = st.text_input(label='Enter Mole Fractions decomp',
										 placeholder="For example: 0.2-0.3-0.5",
										 label_visibility='visible')
			if mol_user_inp != '':
				mol_user_inp_list = np.array(mol_user_inp.split('-')).astype(float)
			else:
				mol_user_inp_list = []

			print(mol_user_inp_list)
			if len(mol_user_inp_list) < len(ele_list_user_inp):
				st.write(":red[Invalid Input]")
			elif len(mol_user_inp_list) <= 1:
				st.write(":red[No Input]")
			elif sum(mol_user_inp_list) != 1:
				st.write(":red[mole fractions must always sum upto 1")
			add_temp_user_inp = st.text_input(label='Enter Temperature',
											placeholder="For example: 1000",
											label_visibility='visible')


			st.write(f"Decomposition products at {add_temp_user_inp} K")
			st.write(pD.find_decomp_products(composition=ele_list_user_inp,
											 mol_ratio=mol_user_inp_list,
											 temperature=add_temp_user_inp,
											 lattice='min')[0])


