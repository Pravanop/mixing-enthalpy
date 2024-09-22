

Plot Polar Plots
===========================

.. code-block :: python


    from far_heaa.io.metadata_handler import MetadataHandler
    from far_heaa.visualization_toolkit.polar_visualizations import PolarVisualizations


    mH = MetadataHandler()
    meta_data = mH.get_metadata

    pD = PolarVisualizations(
        composition=["Cr", "Ta", "W"], lattice="min", meta_data=meta_data, save_flag=True
    )

    _, _ = pD.plot_total()

    _, _ = pD.plot_subset(N_ind=1, transmute_indices=[1, 2])


