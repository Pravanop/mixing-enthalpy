from calcEnthalpy_package.visualization_toolkit.miscellaneous_visualizations import MiscellaneousVisualizations


def mock_meta_data():
	meta_data =  {
		'grid_size': 15,
		'folder_path': '/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/app_data',
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

bD = MiscellaneousVisualizations(
	lattice='min',
	meta_data=meta_data,
	save_flag=True)

ax, fig = bD.plot_convex_hull(composition=['Cr', 'W', 'Ti'],
							  temperature=1500)
