import pytest
import numpy as np
from matplotlib import pyplot as plt
from calcEnthalpy_package.math_operations.thermo_calculations import ThermoMaths
from calcEnthalpy_package.phase_diagram.grid_iterators import GridIterator
from calcEnthalpy_package.visualization_toolkit.binary_visualizations import binaryVizualization


# Mocked data for meta_data
@pytest.fixture
def mock_meta_data():
	meta_data =  {
		'grid_size': 10,
		'folder_path': '/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/app_data',
		'file_name': {
		'unbiased': 'all_lattices_binaries',
		'biased': 'bokas_omegas_processed'
	},
		'end_member': 'bokas_end_members_dict',
		'api_key': 'dummy_api_key',
		'flags': {
			'correction': True,
			'im_flag': True,
			'equi_flag': False
		}
	}
	return meta_data

# Test the initialization of binaryVizualization class
def test_initialization(mock_meta_data):
	composition = ["Fe", "Ni"]
	lattice = "BCC"
	
	bv = binaryVizualization(composition=composition, lattice=lattice, meta_data=mock_meta_data)
	
	# Check if attributes are initialized correctly
	assert bv.composition == composition
	assert bv.lattice == lattice
	assert bv.mol_grid_size == 30
	assert isinstance(bv.tm, ThermoMaths)
	assert isinstance(bv.grid_iterator, GridIterator)


# Test the miscibility temperature function
def test_find_misc_temperatures(mock_meta_data, monkeypatch):
	composition = ["Fe", "Ni"]
	lattice = "BCC"
	
	bv = binaryVizualization(composition=composition, lattice=lattice, meta_data=mock_meta_data)
	
	# Mock the output of grid_iterator.misc_temperature_across_grid
	def mock_misc_temperature_across_grid(composition, mol_grid_size, lattice, phase_flag):
		# Returning dummy mol_grid and misc_temp data
		mol_grid = np.linspace(0, 1, mol_grid_size)
		misc_temp = np.linspace(1000, 1500, mol_grid_size)
		return mol_grid.reshape(-1, 1), misc_temp
	
	monkeypatch.setattr(bv.grid_iterator, "misc_temperature_across_grid", mock_misc_temperature_across_grid)
	
	mol_grid, misc_temp = bv.find_misc_temperatures()
	
	# Validate the results
	assert mol_grid.shape == (30, 1)
	assert len(misc_temp) == 30
	assert np.all(mol_grid >= 0) and np.all(mol_grid <= 1)
	assert np.all(misc_temp >= 1000) and np.all(misc_temp <= 1500)


# Test the plot_misc_temperatures function
def test_plot_misc_temperatures(mock_meta_data, monkeypatch):
	composition = ["Fe", "Ni"]
	lattice = "BCC"
	
	bv = binaryVizualization(composition=composition, lattice=lattice, meta_data=mock_meta_data)
	
	# Mock find_misc_temperatures
	def mock_find_misc_temperatures():
		mol_grid = np.linspace(0, 1, 30).reshape(-1, 1)
		misc_temp = np.linspace(1000, 1500, 30)
		return mol_grid, misc_temp
	
	
	
	monkeypatch.setattr(bv, "find_misc_temperatures", mock_find_misc_temperatures)
	# monkeypatch.setattr(bv.tm, "avg_T_melt", mock_avg_T_melt)
	
	ax, fig = bv.plot_misc_temperatures()
	
	# Check that the figure and axes are created
	assert isinstance(ax, plt.Axes)
	assert isinstance(fig, plt.Figure)
	
	# Check that the plot title is correctly set
	assert ax.get_title() == "Ni-Fe Phase Diagram"
	
	# Ensure the x and y limits are properly set
	# assert ax.get_xlim() == (0, 1)
	# assert ax.get_ylim() == (0, 1800 + 100)
	
	# Check labels
	assert ax.get_xlabel() == 'x'
	assert ax.get_ylabel() == 'T (K)'


# Run the tests
if __name__ == "__main__":
	pytest.main()
