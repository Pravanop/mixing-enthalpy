import streamlit as st
import sys
import os

sys.path.insert(0, os.path.abspath('../..'))

from far_heaa.visualization_toolkit.miscellaneous_visualizations import MiscellaneousVisualizations


def convex_hull(input_list, meta_data):
	if len(input_list) > 3:
		st.write("Convex Hulls are only supported for 3 or fewer elements!")
	T = st.slider('Temperature', 0, 3000, step=200)
	mV = MiscellaneousVisualizations(meta_data=meta_data, save_flag=False)
	ax, fig = mV.plot_convex_hull(input_list, T)
	return ax, fig
