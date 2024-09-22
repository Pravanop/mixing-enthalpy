from far_heaa.visualization_toolkit.heatmap_visualizations import MatrixHeatmap


from far_heaa.io.metadata_handler import MetadataHandler

mH = MetadataHandler()
meta_data = mH.get_metadata

mH = MatrixHeatmap(
    composition=["W", "Cr", "V"],
    add_ele=["W", "Ta"],
    lattice="min",
    meta_data=meta_data,
    save_flag=True,
    path_type="transmutate",
)

mH.plot_ehull_matrix()
