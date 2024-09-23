import streamlit as st
from input_check import InputCheck, get_metadata
from phase_diagrams import phase_diagram_visualizations
from high_symmetry_paths import high_symmetry_paths
from deposition_pathways import deposition_pathways
from convex_hull import convex_hull


# TODO: Add high symmetry polar plots to high symmetry paths
# TODO: Add polar plots for convex hulls
# TODO: Add miscibility temperatures 
st.title('Fast and Robust High Entropy Alloy Analysis')
# take input in specific manner
input_str = st.text_input(label='main_input',
						  placeholder='For example: Al-Cu-Ni',
						  label_visibility='hidden')

# perform sanity checks on the input
if input_str:
	input_check = InputCheck(input_str)
	invalid_flag = input_check.input_validation()
	input_list = input_check.input_list
else:
	invalid_flag = False
	input_list = []

if not invalid_flag:
	
	meta_data = get_metadata()
	
	st.write("What would you like to do?")
	# select the analysis type
	analysis_type = st.selectbox('Analysis type',
								 ['Phase Diagrams', 'High Symmetry Paths', 'Convex Hulls', 'Deposition Pathways'])
	
	if analysis_type == 'Phase Diagrams':
		try:
			ax, fig = phase_diagram_visualizations(input_list, meta_data, lattice='min')
			st.pyplot(fig)
		except TypeError:
			pass
	elif analysis_type == 'High Symmetry Paths':
		try:
			ax, fig = high_symmetry_paths(input_list, meta_data, lattice='min')
			st.pyplot(fig)
		except TypeError:
			pass
	
	elif analysis_type == 'Convex Hulls':
		st.write("Convex Hulls")
		try:
			ax, fig = convex_hull(input_list, meta_data)
			st.pyplot(fig)
		except TypeError:
			pass
	
	elif analysis_type == 'Deposition Pathways':
		st.write("Deposition Pathways")
		try:
			ax, fig = deposition_pathways(input_list, meta_data, lattice='min')
			st.pyplot(fig)
		except TypeError:
			pass
