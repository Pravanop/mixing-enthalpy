from far_heaa.io.metadata_handler import MetadataHandler
from far_heaa.visualization_toolkit.polar_visualizations import PolarVisualizations


mH = MetadataHandler()
meta_data = mH.get_metadata

pD = PolarVisualizations(
    composition=["Al", "Ru", "Cu", "Ni"], lattice="min", meta_data=meta_data, save_flag=True, type_flag='e_hull'
)

_, _ = pD.plot_total(temperature=3000)

# _, _ = pD.plot_subset(N_ind=1, transmute_indices=[], temperature=3000)
