import streamlit as st
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))

from far_heaa.visualization_toolkit.visualizations import Visualizations


def composition_information(input_list, meta_data, lattice):
    sub_options = st.selectbox(
        "What would you like to do?",
        options=["Miscible Temperature", "Decomposition Products"],
    )

    viz = Visualizations(lattice=lattice, meta_data=meta_data)

    gi = viz.grid_iterator

    mol_inputs = []
    cols = st.columns(len(input_list))
    for i in range(len(input_list)):
        with cols[i]:
            mol_inputs.append(
                st.number_input(
                    f"Mole fraction of {input_list[i]}",
                    min_value=0.0,
                    max_value=1.0,
                    value=0.0,
                )
            )

    if sub_options == "Miscible Temperature":
        sum_check = sum(mol_inputs) == 1.0
        if sum_check:
            with st.spinner("Wait for it..."):
                misc_T = gi.uni_molar_misc_temperature(
                    mol_ratio=mol_inputs,
                    composition=input_list,
                    lattice=lattice,
                )
                if isinstance(misc_T, tuple):
                    st.subheader(f"Not Miscible. Alloy melts at {int(misc_T[0])} K.")
                else:
                    st.subheader(f"Miscible temperature: {int(misc_T)} K")
        else:
            st.write("Sum of mole fractions must equal 1.0")

    if sub_options == "Decomposition Products":

        sum_check = sum(mol_inputs) == 1.0

        if sum_check:
            temperature = st.number_input(
                "Temperature (K)", min_value=0, max_value=3800, value=0
            )
            with st.spinner("Wait for it..."):
                decomp, e_hull = gi.find_decomp_products(
                    mol_ratio=mol_inputs,
                    composition=input_list,
                    lattice=lattice,
                    temperature=temperature,
                )
            st.dataframe(decomp)
            # st.write(decomp)
            st.write(f"Energy above convex hull: {e_hull * 1000} meV/atom")
        else:
            st.write("Sum of mole fractions must equal 1.0")
