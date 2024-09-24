import streamlit as st
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))

from far_heaa.visualization_toolkit.binary_visualizations import binaryVizualization
from far_heaa.visualization_toolkit.ternary_visualizations import TernaryVisualization
from far_heaa.visualization_toolkit.polar_visualizations import PolarVisualizations

if "clicked" not in st.session_state:
    st.session_state.clicked = False


def click_button():
    st.session_state.clicked = True


def phase_diagram_visualizations(input_list, meta_data, lattice):
    n_alloy = len(input_list)
    if n_alloy == 2:
        st.write("Binary Phase Diagrams")
        bD = binaryVizualization(
            composition=input_list,
            lattice=lattice,
            meta_data=meta_data,
            save_flag=False,
        )
        with st.spinner("Wait for it..."):
            ax, fig = bD.plot_misc_temperatures()
        return ax, fig
    elif n_alloy == 3:
        st.write("Ternary Phase Diagrams")
        plot_options = st.radio(
            "Select the type of Phase Diagram",
            ["Miscible T contours", "Isotherms"],
            index=None,
        )
        if plot_options:
            tD = TernaryVisualization(
                composition=input_list,
                lattice=lattice,
                meta_data=meta_data,
                save_flag=False,
                contour_flag=True,
            )
            if plot_options == "Miscible T contours":
                with st.spinner("Wait for it..."):
                    ax, fig = tD.plot_misc_temperatures()
            elif plot_options == "Isotherms":
                T = st.slider(
                    "Select the temperature",
                    min_value=0,
                    max_value=3000,
                    value=300,
                    step=300,
                )
                with st.spinner("Wait for it..."):
                    ax, fig = tD.plot_isotherm(T)

            return ax, fig

    elif n_alloy > 3:
        with st.spinner("Wait for it..."):
            pV = PolarVisualizations(
                composition=input_list,
                lattice=lattice,
                meta_data=meta_data,
                save_flag=False,
                type_flag="misc_T",
            )
            with st.spinner("Wait for it..."):
                ax, fig = pV.plot_total()
            return ax, fig
    else:
        st.write("Please provide more elements!")
