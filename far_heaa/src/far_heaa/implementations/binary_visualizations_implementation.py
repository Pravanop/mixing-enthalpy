from far_heaa.visualization_toolkit.binary_visualizations import binaryVizualization
from far_heaa.io.metadata_handler import MetadataHandler

mH = MetadataHandler()
meta_data = mH.get_metadata

bD = binaryVizualization(
    composition=["Ta", "Ti"], lattice="min", meta_data=meta_data, save_flag=True
)

ax, fig = bD.plot_misc_temperatures()
