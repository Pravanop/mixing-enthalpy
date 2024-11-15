from far_heaa.visualization_toolkit.heatmap_visualizations import MatrixHeatmap
from far_heaa.io.metadata_handler import MetadataHandler

mH = MetadataHandler()
meta_data = mH.get_metadata
# meta_data["flags"]["equi_flag"] = True
# meta_data["flags"]["correction"] = False
mh = MatrixHeatmap(
    composition=["Nb", "V", "Zr"],
    add_ele="Ti",
    lattice="min",
    meta_data=meta_data,
    save_flag=True,
    path_type="add",
)

mh.plot_ehull_matrix()
