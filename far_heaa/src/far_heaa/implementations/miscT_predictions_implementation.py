from far_heaa.high_throughput.equi_miscT_predictions import EquiMiscTPredictions
from far_heaa.io.metadata_handler import MetadataHandler

mH = MetadataHandler()
meta_data = mH.get_metadata
# meta_data["file_name"]["biased"] = "bokas_omegas_processed_subregular"
# meta_data["flags"]["equi_flag"] = True
bD = EquiMiscTPredictions(
    element_list=['Ta','V','W','Ti', 'Cr', 'Fe', 'Hf', 'Mn'],
    lattice="BCC",
    meta_data=meta_data,
    save_flag=True,
    generate_plot=True,
    dim=5,
)

ax, fig = bD.plot_predictions()
