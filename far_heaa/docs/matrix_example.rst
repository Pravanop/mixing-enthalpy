

Plot Ehull Diagrams of High Symmetry Paths
===========================================

.. code-block :: python


    from far_heaa.visualization_toolkit.heatmap_visualizations import MatrixHeatmap
    from far_heaa.io.metadata_handler import MetadataHandler

    mH = MetadataHandler()
    meta_data = mH.get_metadata
    # High symmetry paths between equimolar compositions
    mH = MatrixHeatmap(composition=['W', 'Cr', 'V'],
                       add_ele=['W','Ta'],  #first index is starting element and second element is the one it transmutates to
                       lattice='min',
                       meta_data=meta_data,
                       save_flag=True,
                       path_type='transmutate') #if you want a transmutation path

    mH = MatrixHeatmap(composition=['W', 'Cr', 'V'],
                       add_ele=['Ti'], #adds Ti to the ternary above with
                       lattice='min',
                       meta_data=meta_data,
                       save_flag=True,
                       path_type='add') #if you want to add a path

    mH.plot_ehull_matrix()


