

Plot Ternary Phase Features
===========================

.. code-block :: python


    from far_heaa.visualization_toolkit.ternary_visualizations import TernaryVizualization
    from far_heaa.io.metadata_handler import MetadataHandler

    mH = MetadataHandler()
    meta_data = mH.get_metadata

    tD = TernaryVizualization(
        composition=['Cr', 'V', 'Ta'],
        lattice='min',
        meta_data=meta_data,
        save_flag=True,
        contour_flag=False)


    # Plot Miscible Temperatures
    _, _ = tD.plot_misc_temperatures()
    # Plot Isotherm at one temperature
    _, _ = tD.plot_isotherm(1300)

    # Plot contours, scatters for miscible temperatures and defined isotherms
    tD.plot_ternary_visualizations(T_min=300, T_max=3000, T_gradation=500)


