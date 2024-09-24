import streamlit as st
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))

from far_heaa.high_throughput.equi_miscT_predictions import EquiMiscTPredictions


def high_throughput_screening(input_list, meta_data, lattice):
    dim_input = st.number_input(
        "Dimensionality of the system", min_value=2, max_value=10, value=4
    )
    if dim_input:
        len_check = len(input_list) <= dim_input
        if len_check:
            st.write(
                "The dimensionality of the predictions should be lesser than the input."
            )
        else:
            with st.spinner("Wait for it..."):
                equi_mT = EquiMiscTPredictions(
                    dim=dim_input,
                    element_list=input_list,
                    lattice=lattice,
                    meta_data=meta_data,
                )
                fig, ax = equi_mT.plot_predictions()
            return fig, ax

    return None
