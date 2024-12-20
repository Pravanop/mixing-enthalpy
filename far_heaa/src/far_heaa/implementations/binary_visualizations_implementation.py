from far_heaa.visualization_toolkit.binary_visualizations import binaryVizualization
from far_heaa.io.metadata_handler import MetadataHandler

mH = MetadataHandler()
mH.update_metadata(key="flags", value={"im_flag": False, "correction": True, "equi_flag": True})
meta_data = mH.get_metadata
bD = binaryVizualization(
    composition=["Cr", "W"], lattice="min", meta_data=meta_data, save_flag=True
)

ax, fig = bD.plot_misc_temperatures()
