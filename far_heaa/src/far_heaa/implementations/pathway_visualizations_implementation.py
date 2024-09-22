from far_heaa.visualization_toolkit.pathway_visualizations import PathwayVisualizations
from far_heaa.io.metadata_handler import MetadataHandler

mH = MetadataHandler()
meta_data = mH.get_metadata

rP = PathwayVisualizations(
    meta_data=meta_data,
    lattice="BCC",
    composition=["Cr", "V", "Ti", "Ta"],
    save_flag=True,
)

_, _ = rP.plot_rP()
