from calcEnthalpy_package.visualization_toolkit.polar_visualizations import PolarVisualizations


def mock_meta_data():
	meta_data =  {
		'grid_size': 15,
		'folder_path': '../database',
		'file_name': {
		'unbiased': 'all_lattices_binaries',
		'biased': 'bokas_omegas_processed'
	},
		'end_member': 'bokas_end_members_dict',
		'api_key': 'u1TjwfwfTnpF8IolXF9PBY9RT9YauL84',
		'flags': {
			'correction': True,
			'im_flag': True,
			'equi_flag': False
		}
	}
	return meta_data


meta_data = mock_meta_data()

pD = PolarVisualizations(
	composition=['Cr', 'Ta', 'W', 'V','Ti'],
	lattice='min',
	meta_data=meta_data,
	save_flag=True)

_, _ = pD.total_plot()

_, _ = pD.plot_N(N_ind=1,
				 transmute_indices=[1, 2])


