from calcEnthalpy_package.visualization_toolkit.pathway_visualizations import PathwayVisualizations

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

rP = PathwayVisualizations(
	meta_data=meta_data,
	lattice='BCC',
	composition=['Cr', 'V', 'Ti', 'Ta','W'],
	save_flag=True
)

_, _ = rP.plot_rP()