from far_heaa.visualization_toolkit.miscellaneous_visualizations import (
    MiscellaneousVisualizations,
)
from far_heaa.io.metadata_handler import MetadataHandler

mH = MetadataHandler()
meta_data = mH.get_metadata

bD = MiscellaneousVisualizations(meta_data=meta_data, save_flag=True)

ax, fig = bD.plot_convex_hull(composition=["Cr", "W", "Ti"], temperature=1500)
