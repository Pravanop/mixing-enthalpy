
from far_heaa.io.metadata_handler import MetadataHandler
from far_heaa.visualization_toolkit.polar_visualizations import PolarVisualizations


mH = MetadataHandler()
meta_data = mH.get_metadata
meta_data['flags']['equi_flag'] = True
# meta_data['flags']['correction'] = False
# meta_data['flags']['im_flag'] = False
meta_data['file_name']['biased'] = 'bokas_old_omegas_processed'
pD = PolarVisualizations(
    # composition=["Cr", "Co", "Mn", "Fe", "Ni"],
    # composition=["Nb", "Mo", "Hf", "Ti", "Zr"],
    composition=["Nb", "Mo", "Hf"],
    # composition=["A", "B", "D", "C"],
    # composition=["Al", "Co", "Cu", "Fe", "Ni"],
    lattice="min",
    meta_data=meta_data,
    save_flag=True,
    type_flag="e_hull",
)

_, _ = pD.plot_total(temperature=563)

#
# # # # #
# _, _ = pD.plot_subset(N_ind=2,
#                       transmute_indices=[])
# # # #
# _, _ = pD.plot_subset(N_ind=1,
#                       transmute_indices=[1, 2])

