from far_heaa.visualization_toolkit.binary_visualizations import binaryVizualization
from far_heaa.io.metadata_handler import MetadataHandler
from far_heaa.visualization_toolkit.ternary_visualizations import TernaryVisualization
from far_heaa.visualization_toolkit.heatmap_visualizations import MatrixHeatmap
from far_heaa.visualization_toolkit.miscellaneous_visualizations import (
    MiscellaneousVisualizations,
)
from far_heaa.visualization_toolkit.pathway_visualizations import PathwayVisualizations
from far_heaa.visualization_toolkit.polar_visualizations import PolarVisualizations

# Load metadata
mH = MetadataHandler()
meta_data = mH.get_metadata

# Binary visualization
bD = binaryVizualization(
    composition=["Ta", "Ti"], lattice="min", meta_data=meta_data, save_flag=True
)

_, _ = bD.plot_misc_temperatures()

# # Heatmap visualization
mH = MatrixHeatmap(
    composition=["W", "Cr", "V"],
    add_ele=["W", "Ta"],
    lattice="min",
    meta_data=meta_data,
    save_flag=True,
    path_type="transmutate",
)

mH.plot_ehull_matrix()
#
#
# # Miscellaneous visualization
mV = MiscellaneousVisualizations(meta_data=meta_data, save_flag=True)

_, _ = mV.plot_convex_hull(composition=["Cr", "W", "Ti"], temperature=1500)
#
#
# # Pathway visualization
rP = PathwayVisualizations(
    meta_data=meta_data,
    lattice="BCC",
    composition=["Cr", "V", "Ti", "Ta"],
    save_flag=True,
)

_, _ = rP.plot_rP()
#
#
# # Polar visualization
pD = PolarVisualizations(
    composition=["Cr", "Ta", "W", "V"],
    lattice="min",
    meta_data=meta_data,
    save_flag=True,
)

_, _ = pD.plot_total()

_, _ = pD.plot_subset(N_ind=1, transmute_indices=[1, 2])
#
#
# # Ternary visualization
tD = TernaryVisualization(
    composition=["Cr", "V", "Ta"],
    lattice="min",
    meta_data=meta_data,
    save_flag=True,
    contour_flag=False,
)

# Plot everything
tD.plot_ternary_visualizations(T_min=300, T_max=3000, T_gradation=500)
