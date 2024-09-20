import os

import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm

from calcEnthalpy_package.grids_and_combinations.combination_generation import MultinaryCombinations
from calcEnthalpy_package.io.json_handler import JSONHandler
from calcEnthalpy_package.math_operations.thermo_calculations import ThermoMaths
from calcEnthalpy_package.phase_diagram.grid_iterators import GridIterator


class EquiMiscTPredictions:
	
	def __init__(self,
				 dim,
				 element_list,
				 lattice,
				 meta_data
				 ):
		
		self.lattice = lattice
		self.mol_grid_size = 40
		
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
		
		self.dim = dim
		self.element_list = element_list
		
		compositions = MultinaryCombinations.create_multinary(element_list=element_list,
															  no_comb=[self.dim])
		
		self.compositions = list(compositions.values())[0]
		self.tm = ThermoMaths()
	
	def make_predictions(self):
		
		avg_tm = []
		temp_list = []
		alloys = []
		for composition in tqdm(self.compositions, desc='Calculating Misc_T'):
			composition = composition.split('-')
			n_alloy = len(composition)
			
			mol = [1 /n_alloy] * n_alloy
			misc_T = self.grid_iterator.uni_molar_misc_temperature(
				mol_ratio=mol,
				composition=composition,
				lattice=self.lattice,
				phase_flag=True,
				batch_tag=False
			)
			if isinstance(misc_T, float):
				temp_list.append(misc_T)
			else:
				temp_list.append(misc_T[0])
			alloys.append('-'.join(composition))
			avg_tm.append(self.tm.avg_T_melt(composition=composition,
											 mol_ratio=mol))
		
		df = pd.DataFrame([alloys, temp_list, avg_tm])
		df = df.T
		df.columns = ["Alloys", 'Misc_Temp', 'Avg_Tm']
		
		if not os.path.exists("../output_data/predictions/"):
			os.mkdir("../output_data/predictions")
		df.to_csv(f"../output_data/predictions/misc_T_{self.dim}_{self.lattice}.csv", index=False)
		
		return df
	
	def plot_predictions(self,
						 Tm_constraint,
						 alloy_constraint
							):
		file_path = "../output_data/predictions/misc_T_{self.dim}_{self.lattice}.csv"
		
		if not os.path.exists(file_path):
			df = self.make_predictions()
		else:
			df = pd.read_csv(f"./misc_T_{self.dim}_equi_{self.lattice}.csv")
		
		if Tm_constraint is not None:
			df_proc = df.loc[df['Temp'] < Tm_constraint * df['Avg_Tm']]
		if alloy_constraint is not None:
			df_proc = df_proc[df['Alloys'].str.contains(alloy_constraint)]
		
		if Tm_constraint is None and alloy_constraint is None:
			df_proc =  df
			
		df_proc = df_proc.sort_values(['Avg_Tm'])
		fig, ax = plt.subplots(figsize=(18, 8))
		ax.bar(x=df_proc['Alloys'], height=df_proc['Avg_Tm'], color='#4393C3', edgecolor='black', width=0.5,
			   align='center', zorder=2)
		ax.bar(x=df_proc['Alloys'], height=df_proc['Temp'], color='#D6604D', edgecolor='black', width=0.5,
			   align='center', zorder=2)
		ax.axhline(y=2000, color='black', linestyle='--', linewidth=3.5, alpha=0.8, zorder=0, label='_nolegend_')
		ax.legend(['Average T$_{m}$', 'Miscible Temperature', ])
		# ax.set_xlabel('Alloys')
		ax.set_ylabel('T (K)')
		ax.tick_params(axis='x', labelrotation=90)
		plt.subplots_adjust(bottom=0.32, top=0.99, left=0.1, right=0.99)
		
		if not os.path.exists("../plots/predictions"):
			os.mkdir("../plots/predictions")
			
		plt.savefig(f'../plots/predictions/misc_T_{self.dim}_{self.lattice}_{Tm_constraint}_{alloy_constraint}.png')
