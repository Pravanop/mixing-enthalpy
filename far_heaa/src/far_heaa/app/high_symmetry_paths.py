import streamlit as st
import sys
import os
from input_check import InputCheck

sys.path.insert(0, os.path.abspath("../.."))

from far_heaa.visualization_toolkit.polar_visualizations import PolarVisualizations


def high_symmetry_paths(input_list, meta_data, lattice):

    if len(input_list) <= 3:
        st.write(
            "Why don't you just plot the phase diagram? Or provide a higher order phase space."
        )
        return None

    st.write(
        f"If you pass N_ind=1, N_ind={len(input_list)-1} will also be used. Thus, it is wise to pass the lower index."
    )
    with st.spinner("Wait for it..."):
        pD = PolarVisualizations(
            composition=input_list,
            lattice=lattice,
            meta_data=meta_data,
            save_flag=False,
            type_flag="misc_T",
        )
        N_ind = st.number_input(
            "Index of High Symmetry Path (1 is to unary, 2 to binary and so on)",
            min_value=1,
            max_value=len(input_list) - 1,
            value=1,
        )
        ax, fig = pD.plot_subset(N_ind=N_ind, transmute_indices=[])

    return ax, fig
