import streamlit as st
import sys
import os
sys.path.insert(0, os.path.abspath('../..'))


from far_heaa.visualization_toolkit.binary_visualizations import binaryVizualization
from far_heaa.visualization_toolkit.ternary_visualizations import TernaryVisualization

def phase_diagram_visualizations(input_list, meta_data, lattice):
	
	n_alloy = len(input_list)
	if n_alloy == 2:
		st.write("Binary Phase Diagrams")
		bD = binaryVizualization(
			composition=input_list, lattice=lattice, meta_data=meta_data, save_flag=False
		)
		ax, fig = bD.plot_misc_temperatures()
		return ax, fig
	elif n_alloy == 3:
		st.write("Ternary Phase Diagrams")
		plot_options = st.radio("Select the type of Phase Diagram", ['Miscible T contours', 'Isotherms'])
		st.button("Plot")
		tD = TernaryVisualization(
			composition=input_list, lattice=lattice, meta_data=meta_data, save_flag=False, contour_flag=True
		)
		if plot_options == 'Miscible T contours':
			ax, fig = tD.plot_misc_temperatures()
		elif plot_options == 'Isotherms':
			T = st.slider("Select the temperature", min_value=0, max_value=3000, value=300, step = 300)
			ax, fig = tD.plot_isotherm(T)
			
		return ax, fig
		
	elif n_alloy > 3:
		st.write("Multi-component Phase Diagrams")
	else:
		st.write("Please provide more elements!")
	
	