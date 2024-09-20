import matplotlib.pyplot as plt
from pymatgen.analysis.phase_diagram import PDPlotter
import matplotlib as mpl
from calcEnthalpy_package.io.dir_handler import DirHandler
from calcEnthalpy_package.io.json_handler import JSONHandler
from calcEnthalpy_package.math_operations.thermo_calculations import ThermoMaths
from calcEnthalpy_package.phase_diagram.grid_iterators import GridIterator


class MiscellaneousVisualizations:
	
	def __init__(self,
				 lattice,
				 meta_data,
				 save_flag):
		
		self.lattice = lattice
		self.mol_grid_size = 40
		self.flags = meta_data['flags']
		self.tm = ThermoMaths()
		grid_size = meta_data['grid_size']
		
		if meta_data['flags']['correction']:
			data = JSONHandler.load_json(folder_path=meta_data['folder_path'],
										 file_name=meta_data['file_name']['biased'])
		else:
			data = JSONHandler.load_json(folder_path=meta_data['folder_path'],
										 file_name=meta_data['file_name']['unbiased'])
		
		end_member = JSONHandler.load_json(folder_path=meta_data['folder_path'], file_name=meta_data['end_member'])
		
		self.grid_iterator = GridIterator(grid_size=grid_size,
										  tm=self.tm,
										  data=data,
										  end_member=end_member,
										  api_key=meta_data['api_key'],
										  flags=meta_data['flags']
										  )
		self.save_flag = save_flag
	
	def plot_convex_hull(self, composition, temperature):
		
		conv_hull = self.grid_iterator.convex_hull.make_convex_hull(temperature=temperature,
																	composition=composition,
																	batch_tag=False)
		
		fig, ax = plt.subplots(figsize=(3.2, 3.2))
		mpl.rcParams.update({'font.size': 10})
		PDPlotter(phasediagram=conv_hull,
				  show_unstable=0.02,
				  backend='matplotlib').get_plot(
			label_unstable=False,
			label_stable=False,
			process_attributes=True,
			ax=ax,
			fill=False
		)
		if len(composition) == 2:
			ax.set_xlim([0, 1])
			ax.set_xlabel('x', fontsize=10)
			ax.set_ylabel('E hull (eV/atom)', fontsize=10)
			plt.subplots_adjust(bottom=0.12, left=0.15, right=0.95, top=0.85)
			ax.set_title(f'{"".join(sorted(composition, reverse = True))}')
			
		if len(composition) == 3:
			plt.subplots_adjust(top=0.95)
		
		vert = [(0.5, 0.9), (-0.05, 0), (1.02, 0)]
		composition.reverse()
		for idx, i in enumerate(composition):
			ax.text(s = i, x = vert[idx][0], y = vert[idx][1], fontsize = 12)
		if self.save_flag:
			updated_folder_path = DirHandler.mkdir_recursrive(
				folders=['convex_hull', f'{len(composition)}', f'{"".join(sorted(composition))}'],
				folder_path='../plots')
			
			if self.flags['equi_flag']:
				fig.savefig(fname=f"{updated_folder_path}/convex_hull_{temperature}_equi.png", dpi=100)
			else:
				fig.savefig(fname=f"{updated_folder_path}/convex_hull_{temperature}.png", dpi=100)
		
		return ax, fig
