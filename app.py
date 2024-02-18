import streamlit as st
from UI.search_intermetallics import search
from pymatgen.analysis.phase_diagram import PDPlotter
import pandas as pd

def convert_df(df) :
	# IMPORTANT: Cache the conversion to prevent computation on every rerun
	return df.to_csv().encode('utf-8')

st.title('High Entropy Alloy Analysis from First-Principles (HEA-AFP)')

user_inp = st.text_input(label = "A list of elements in form Al-Cu-Ni")
ele_list_user_inp = user_inp.split('-')
invalid = True
if user_inp:
	invalid = False
	try :
		search_inst = search(inp = ele_list_user_inp)
	except ValueError as e :
		invalid = True
		st.write(str(e))
else:
	st.write("Give an input to use this database!")

if not invalid :
	col1 , col2 , col3 , col4 = st.columns([1 , 1 , 1 , 1])
	with col2 :
		find_enthalpy = st.button("Find HEA Enthalpy")
	with col3 :
		find_intermetallic = st.button("Find intermetallics")
	
	if find_enthalpy :
		answer = search_inst.search_equimolar
		st.write(answer)
		st.download_button("Download" , answer , file_name = f"{user_inp}_enthalpy.txt")
	
	if find_intermetallic :
		answer = search_inst.search_all_intermetallics
		st.write(answer)
		if not isinstance(answer , str) :
			df_answer = convert_df(answer)
			st.download_button("Download as csv" , df_answer , file_name = f"{user_inp}.csv")
	
	combinations = st.checkbox("Find N-nary Combinations from the provided composition")
	if combinations :
		st.write("You can check multiple combinations too...")
		opts = ["binary" , "ternary" , "quaternary" , "quinary"]
		known_variables = {symbol : st.checkbox(symbol) for symbol in opts}
		comb_done = st.button("Done")
		
		if comb_done :
			no_combinations = [key if value is True else None for key , value in known_variables.items()]
			
			answer = search_inst.search_combinations(no_combinations)
			st.write(answer)
	
	st.write("===========")
	phase_diagram = st.checkbox("Find Phase Diagram")
	if phase_diagram :
		T = st.slider(label = "Temperature (K)" , min_value = 0 , max_value = 2500)
		answer = search_inst.get_phase_diagram(int(T))
		st.write(f"Temperature: {T} K")
		try:
			st.write(PDPlotter(answer , show_unstable = True).get_plot())
		except ValueError as e:
			st.write("Cannot Visualize Phase Diagram! ",str(e))
			st.write("But here are the stable entries at this temperature")
			error_answer = [[i.as_dict()['name'],i.as_dict()['energy']] for i in list(answer.stable_entries)]
			st.write(pd.DataFrame(error_answer, columns = ['Formula', 'Energy']).sort_values(by = 'Energy',
			                                                                                 ignore_index = True))
			]
			