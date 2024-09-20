from calcEnthalpy_package.visualization_toolkit.ternary_visualizations import TernaryVizualization


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

tD = TernaryVizualization(
	composition=['Ti', 'Cr', 'V'],
	lattice='min',
	meta_data=meta_data,
save_flag=True,
contour_flag=False)


"""Two options"""

# Plot Miscible Temperatures
# _, _ = tD.plot_misc_temperatures()
# Plot Isotherms
_, _ = tD.plot_isotherm(1273)
