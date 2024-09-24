import streamlit as st
import sys
import os
from input_check import InputCheck

sys.path.insert(0, os.path.abspath("../.."))

from far_heaa.visualization_toolkit.heatmap_visualizations import MatrixHeatmap


def alloying_paths(input_list, meta_data, lattice):
    st.write("High Symmetry Paths")
    path_type = st.radio("Select the type of path", ["add", "transmutate"], index=False)
    if path_type == "add":
        st.write("Add")
        add_input = st.text_input("Add input", placeholder="For example: Al")
        if add_input:
            input_check = InputCheck(add_input)
            subset_all_flag, subset_database_flag = (
                input_check.subset_all(),
                input_check.subset_database(),
            )

            if not subset_all_flag:
                st.write(
                    f"{', '.join(list(input_check.input_set.difference(input_check.all_element_set)))} not a valid element!"
                )
                return None
            if not subset_database_flag:
                st.write(
                    f"{', '.join(list(input_check.input_set.difference(input_check.database_element_set)))} not a valid element!"
                )
                return None
            if add_input in input_list:
                st.write(
                    f"{add_input} in the original input. Non-stoichiometric composition!"
                )
                return None
            else:
                with st.spinner("Wait for it..."):
                    mH = MatrixHeatmap(
                        meta_data=meta_data,
                        lattice=lattice,
                        add_ele=[add_input],
                        composition=input_list,
                        save_flag=False,
                        path_type=path_type,
                    )
                    ax, fig = mH.plot_ehull_matrix()
                return ax, fig

    elif path_type == "transmutate":
        st.write("Transmutate")
        transmutate_input = st.text_input(
            "Transmutate input", placeholder="For ex: Ta-W"
        )
        if transmutate_input:
            input_check = InputCheck(transmutate_input)
            transmutate_list = input_check.input_list
            subset_all_flag, subset_database_flag = (
                input_check.subset_all(),
                input_check.subset_database(),
            )
            if not subset_all_flag:
                st.write(
                    f"{', '.join(list(input_check.input_set.difference(input_check.all_element_set)))} not a valid element!"
                )
                return None
            if not subset_database_flag:
                st.write(
                    f"{', '.join(list(input_check.input_set.difference(input_check.database_element_set)))} not a valid element!"
                )
                return None
            if transmutate_list[0] not in input_list:
                st.write(f"{transmutate_list[1]} not in the original input")
                return None
            else:
                with st.spinner("Wait for it..."):
                    mH = MatrixHeatmap(
                        meta_data=meta_data,
                        lattice=lattice,
                        add_ele=transmutate_list,
                        composition=input_list,
                        save_flag=False,
                        path_type=path_type,
                    )
                    ax, fig = mH.plot_ehull_matrix()
                return ax, fig

        else:
            return None
