from far_heaa.high_throughput.equi_miscT_predictions import EquiMiscTPredictions
from far_heaa.io.metadata_handler import MetadataHandler

mH = MetadataHandler()
meta_data = mH.get_metadata

bD = EquiMiscTPredictions(
    element_list=["Cr", "Ti", "Ta", "W"],
    lattice="min",
    meta_data=meta_data,
    save_flag=False,
    generate_plot=True,
    dim=3,
)

ax, fig = bD.plot_predictions()
