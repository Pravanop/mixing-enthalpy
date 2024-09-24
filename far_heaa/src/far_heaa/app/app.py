import streamlit as st
from input_check import InputCheck, get_metadata
from phase_diagrams import phase_diagram_visualizations
from alloying_paths import alloying_paths
from deposition_pathways import deposition_pathways
from convex_hull import convex_hull
from high_symmetry_paths import high_symmetry_paths
from high_throughput_screening import high_throughput_screening
from composition import composition_information

# TODO: Add polar plots for convex hulls
# TODO: Add flag management

st.title("Fast and Robust High Entropy Alloy Analysis")
# take input in specific manner
input_str = st.text_input(
    label="main_input", placeholder="Al-Cu-Ni", label_visibility="hidden"
)


def flag_setter(mH, flags):
    mH.update_metadata(key="flags", value=flags)

    return mH.get_metadata


# perform sanity checks on the input
if input_str:
    input_check = InputCheck(input_str)
    invalid_flag = input_check.input_validation()
    input_list = input_check.input_list
else:
    invalid_flag = True
    input_list = []

if not invalid_flag:

    meta_data, mH = get_metadata()

    st.subheader("Lattice")
    lattice = st.radio(
        "Lattice",
        label_visibility="hidden",
        options=["min", "FCC", "BCC", "HCP"],
        index=None,
        horizontal=True,
    )
    st.subheader("Flags")
    flagoptions = st.radio(
        "Flags",
        label_visibility="hidden",
        options=["Without IM", "T Correction off", "Equimolar grid", "No Changes"],
        index=3,
        horizontal=True,
    )
    flags = {"im_flag": True, "correction": True, "equi_flag": False}
    if flagoptions == "Without IM":
        flags["im_flag"] = False
    elif flagoptions == "T Correction off":
        flags["correction"] = False
    elif flagoptions == "Equimolar grid":
        flags["equi_flag"] = True
    elif flagoptions == "No Changes":
        pass

    meta_data = flag_setter(mH, flags)
    st.header("Visualizations")
    analysis_type = st.selectbox(
        "Analysis type",
        label_visibility="hidden",
        options=[
            "Specific Composition Information",
            "Phase Diagrams",
            "High Symmetry Paths",
            "Alloying Paths",
            "Convex Hulls",
            "Deposition Pathways",
            "High Throughput Screening",
        ],
    )
    if lattice:
        if analysis_type == "Specific Composition Information":
            composition_information(input_list, meta_data, lattice=lattice)

        elif analysis_type == "Phase Diagrams":
            try:
                ax, fig = phase_diagram_visualizations(
                    input_list, meta_data, lattice=lattice
                )
                st.pyplot(fig)
            except TypeError:
                pass
        elif analysis_type == "High Symmetry Paths":
            try:
                ax, fig = high_symmetry_paths(input_list, meta_data, lattice=lattice)
                st.pyplot(fig)
            except TypeError:
                pass
        elif analysis_type == "Alloying Paths":
            try:
                ax, fig = alloying_paths(input_list, meta_data, lattice=lattice)
                st.pyplot(fig)
            except TypeError:
                pass

        elif analysis_type == "Convex Hulls":
            try:
                ax, fig = convex_hull(input_list, meta_data, lattice=lattice)
                st.pyplot(fig)
            except TypeError:
                pass

        elif analysis_type == "Deposition Pathways":
            try:
                ax, fig = deposition_pathways(input_list, meta_data, lattice=lattice)
                st.pyplot(fig)
            except TypeError:
                pass

        elif analysis_type == "High Throughput Screening":
            try:
                fig, ax = high_throughput_screening(
                    input_list, meta_data, lattice=lattice
                )
                st.pyplot(fig)
            except TypeError:
                pass
