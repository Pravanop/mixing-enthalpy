from far_heaa.visualization_toolkit.miscellaneous_visualizations import (
    MiscellaneousVisualizations,
)
from far_heaa.io.metadata_handler import MetadataHandler
import matplotlib.pyplot as plt
mH = MetadataHandler()
meta_data = mH.get_metadata

bD = MiscellaneousVisualizations(meta_data=meta_data, save_flag=True)

# ax, fig = bD.plot_convex_hull(composition=["Cr", "W", "Ti"], temperature=1500)
fig = bD.heatmap(element_list=['Cr', 'V', 'W', 'Ti', 'Ta', 'Fe', 'Mo', 'Nb', 'Zr', 'Hf'], lattice = 'BCC')
plt.show()

