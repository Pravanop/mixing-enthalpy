
from far_heaa.io.metadata_handler import MetadataHandler
from far_heaa.visualization_toolkit.polar_visualizations import PolarVisualizations


mH = MetadataHandler()
meta_data = mH.get_metadata
meta_data['flags']['equi_flag'] = True
pD = PolarVisualizations(
    # composition=["Hf", "Mo", "Nb", "Ti", "Zr"],
    composition=["A", "B", "D", "C", "E"],
    lattice="BCC",
    meta_data=meta_data,
    save_flag=True,
    type_flag="gibbs",
)

_, _ = pD.plot_total()


# # #
# _, _ = pD.plot_subset(N_ind=2,
#                       transmute_indices=[])
#
_, _ = pD.plot_subset(N_ind=1,
                      transmute_indices=[1, 2])

