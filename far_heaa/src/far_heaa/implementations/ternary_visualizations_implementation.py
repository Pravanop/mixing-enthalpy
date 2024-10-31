import time

from far_heaa.visualization_toolkit.ternary_visualizations import TernaryVisualization
from far_heaa.io.metadata_handler import MetadataHandler

mH = MetadataHandler()
meta_data = mH.get_metadata
# meta_data["flags"]["equi_flag"] = True
tD = TernaryVisualization(
    composition=["Nb", "Zr", "V"],
    lattice="min",
    meta_data=meta_data,
    save_flag=True,
    contour_flag=False,
)

"""Two options"""

start = time.time()
# Plot Miscible Temperatures
_, _ = tD.plot_misc_temperatures()
stop = time.time()
print(stop - start)
# Plot Isotherms
# _, _ = tD.plot_isotherm(1300)

# Plot everything
# tD.plot_ternary_visualizations(T_min=300, T_max=3000, T_gradation=500)
