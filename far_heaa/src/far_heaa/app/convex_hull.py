import streamlit as st
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))

from far_heaa.visualization_toolkit.miscellaneous_visualizations import (
    MiscellaneousVisualizations,
)
from far_heaa.visualization_toolkit.polar_visualizations import PolarVisualizations


def convex_hull(input_list, meta_data, lattice):
    print(5)
    if len(input_list) > 3:
        T = st.slider("Temperature", 0, 3000, step=200, value=200)
        with st.spinner("Wait for it..."):
            pV = PolarVisualizations(
                composition=input_list,
                meta_data=meta_data,
                save_flag=False,
                type_flag="e_hull",
                lattice=lattice,
            )
            ax, fig = pV.plot_total(temperature=T)
        return ax, fig
    elif len(input_list) <= 3:
        print(3)
        with st.spinner("Wait for it..."):
            T = st.slider("Temperature", 0, 3000, step=200)
            mV = MiscellaneousVisualizations(meta_data=meta_data, save_flag=False)
            ax, fig = mV.plot_convex_hull(input_list, T)
        return ax, fig
