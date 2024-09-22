import streamlit as st
from input_check import input_validation, get_metadata
from phase_diagrams import phase_diagram_visualizations
st.title('Fast and Robust High Entropy Alloy Analysis')
# take input in specific manner
input_str = st.text_input(label='main_input',
						  placeholder='For example: Al-Cu-Ni',
						  label_visibility='hidden')

#perform sanity checks on the input
if input_str:
	invalid_flag, input_list = input_validation(input_str)
else:
	invalid_flag = False
	input_list = []
	
if not invalid_flag:
	
	meta_data = get_metadata()
	
	st.write("What would you like to do?")
	# select the analysis type
	analysis_type = st.selectbox('Analysis type', ['Phase Diagrams', 'High Symmetry Paths', 'Convex Hulls', 'Deposition Pathways'])
	
	if analysis_type == 'Phase Diagrams':
		ax, fig = phase_diagram_visualizations(input_list, meta_data, lattice='min')
		st.pyplot(fig)
	elif analysis_type == 'High Symmetry Paths':
		st.write("High Symmetry Paths")
	
	elif analysis_type == 'Convex Hulls':
		st.write("Convex Hulls")
	
	elif analysis_type == 'Deposition Pathways':
		st.write("Deposition Pathways")