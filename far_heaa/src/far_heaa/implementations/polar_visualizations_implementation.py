from far_heaa.io.metadata_handler import MetadataHandler
from far_heaa.visualization_toolkit.polar_visualizations import PolarVisualizations


mH = MetadataHandler()
meta_data = mH.get_metadata

pD = PolarVisualizations(
    composition=["Cr", "W", "Ta", "V", "Ti"], lattice="min", meta_data=meta_data, save_flag=True, type_flag='misc_T'
)

_, _ = pD.plot_total()

_, _ = pD.plot_subset(N_ind=1, transmute_indices=[1, 2])
