import streamlit as st
from pymatgen.core import Element

from UI.search import search_composition , search_any
from pymatgen.analysis.phase_diagram import PDPlotter
import pandas as pd

def convert_df(df) :
	# IMPORTANT: Cache the conversion to prevent computation on every rerun
	return df.to_csv().encode('utf-8')

def clicked(button) :
	st.session_state.clicked[button] = True

if 'clicked' not in st.session_state :
	st.session_state.clicked = {1 : False , 2 : False, 3: False}
st.set_page_config(layout = "wide")
st.title('High Entropy Alloy Analysis from First-Principles (HEA-AFP)')
check = True
user_inp = st.text_input(label = 'Input',placeholder = "Provide a list of elements in form Al-Cu-Ni")
ele_list_user_inp = user_inp.split('-')
element_list_path = ("/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/data/input_data/pravan/element_list_bcc_pravan.txt")

with open(element_list_path) as f :
	element_total_list = f.read()
element_total_list = element_total_list.split(',')
print(ele_list_user_inp)
if ele_list_user_inp != ['']:
	check = set(ele_list_user_inp) <= set(element_total_list)
if not check :
	st.write(
			f":red["
			f":blue[{', '.join(list(set(ele_list_user_inp).difference(set(element_total_list))))}] not in "
			f"database! We are working on increasing our database!]"
			)

check = len(ele_list_user_inp ) == len(set(ele_list_user_inp ))
if check is False:
	st.write(":red[Repeated Element!]")
	
invalid = True
col1, col2 = st.columns([1, 1])
if user_inp and check:
	
	with col1:
		find_any = st.button('Find Compositions containing these elements' , on_click = clicked , args = [1])
	with col2:
		find_comp = st.button('Find Data Specific to this Composition' , on_click = clicked , args = [2])
	
	with col2:
		if st.session_state.clicked[2] :
			print(st.session_state.clicked[2])
			invalid = False
			try :
				search_inst = search_composition(inp = ele_list_user_inp,)
				search_inst2 = search_composition(inp = ele_list_user_inp,)
			except ValueError as e :
				invalid = True
				st.write(str(e))
	with col1:
		if st.session_state.clicked[1] :
			hide = st.write("What alloy systems are you interested in? (You can check multiple)")
			if True:
				sys_opts = ["binaries" , "ternaries" , "quaternaries" , "quinaries"]
				alloy_sys = {symbol : st.checkbox(symbol) for symbol in sys_opts}
				main_done = st.button("Done", key = "main-done", on_click = clicked, args=[3])
				if st.session_state.clicked[3]:
					find = [key if value is True else None for key , value in alloy_sys.items()]
					while None in find :
						find.remove(None)
					try :
						search_inst_any = search_any(inp = ele_list_user_inp , find = find)
						answer = search_inst_any.search_all_contains
						if isinstance(answer, str):
							st.write(answer)
						else:
							answer = pd.DataFrame(answer).sort_values(
								by = 'Enthalpy' , ignore_index =
								True
								)
							st.write(f"There are {len(answer.loc[answer['Enthalpy'] <= 0])} stable alloys for your selection.")
							unstable = st.checkbox("Would you like Unstable Alloys as well?")
							if unstable:
								st.write(answer)
							else:
								st.write(answer.loc[answer['Enthalpy'] <= 0])
					except ValueError as e :
						st.write(str(e))
elif check:
	st.write(":red[Give an input to use this database!]")
elif user_inp:
	st.write("")

with col2:
	if not invalid :
		
		side_col1 , side_col4 = st.columns([1 , 1])
		with side_col1 :
			find_enthalpy = st.button("Find HEA Enthalpy")
		with side_col4 :
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
			comb_done = st.button("Done", key = 'comb-done')
			
			if comb_done :
				no_combinations = [key if value is True else None for key , value in known_variables.items()]
				
				answer = search_inst.search_combinations(no_combinations)
				st.write(answer)
		
		phase_diagram = st.checkbox("Find Phase Diagram")
		if phase_diagram :
			T = st.slider(label = "Temperature (K)" , min_value = 0 , max_value = 3000)
			answer = search_inst.get_phase_diagram(int(T), flag="offequi")
			# answer2 = search_inst2.get_phase_diagram(int(T), flag="equi")
			st.write(f"Temperature: {T} K")
			try :
				answer_col1, answer_col2 = st.columns([1, 1])
				with answer_col1 :
					st.write(PDPlotter(answer , show_unstable = True, ternary_style='3d').get_plot())
				# with answer_col2 :
				# 	st.write(PDPlotter(answer2 , show_unstable = True).get_plot())
			except ValueError as e :
				st.write("Cannot Visualize Phase Diagram! " , str(e))
				st.write("But here are the stable entries at this temperature")
				error_answer = [[i.as_dict()['name'] , i.as_dict()['energy']] for i in list(answer.stable_entries)]
				st.write(
					pd.DataFrame(error_answer , columns = ['Formula' , 'Energy']).sort_values(
						by = 'Energy' ,
						ignore_index = True
						)
					)
