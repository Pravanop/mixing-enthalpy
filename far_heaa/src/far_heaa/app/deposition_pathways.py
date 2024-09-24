import streamlit as st
import sys
import os
from input_check import InputCheck

sys.path.insert(0, os.path.abspath("../.."))

from far_heaa.visualization_toolkit.pathway_visualizations import PathwayVisualizations


def deposition_pathways(input_list, meta_data, lattice):
    with st.spinner("Wait for it..."):
        rP = PathwayVisualizations(
            meta_data=meta_data,
            lattice=lattice,
            composition=input_list,
            save_flag=False,
        )
        fig, ax = rP.plot_rP()
    return fig, ax
