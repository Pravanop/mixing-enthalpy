import time

from far_heaa.visualization_toolkit.ternary_visualizations import TernaryVisualization
from far_heaa.io.metadata_handler import MetadataHandler

mH = MetadataHandler()
meta_data = mH.get_metadata
tD = TernaryVisualization(
    composition=["Nb", "V", "Zr"],
    lattice="BCC",
    meta_data=meta_data,
    save_flag=True,
    contour_flag=False,
)
print(meta_data)
"""Two options"""

# Plot Miscible Temperatures
plot = tD.plot_misc_temperatures()

# Plot Isotherms
# plot = tD.plot_isotherm(1400)
plot.show()
# Plot everything
# tD.plot_ternary_visualizations(T_min=300, T_max=3000, T_gradation=500)
