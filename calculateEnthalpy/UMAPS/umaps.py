"""
UMAP coordinate generator.

See UMAP documentation here: https://umap-learn.readthedocs.io/en/latest/how_umap_works.html

Author: Trevor Hastings
Contact: trevorhastings@tamu.edu
Creation Date: 2024/01/09
License: MIT License

Description:

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt, matplotlib.colors as colors
import time
import umap.umap_ as umap
import json

from calculateEnthalpy.helper_functions.grid_code import create_mol_grid


def change_directory_to_current_py_file():
	try:
		import os
		os.chdir(os.path.split(os.path.abspath(__file__))[0])
	except Exception as e:
		print(f'\nDirectory change or os module error:\n\n{e}')


def read_data(file_name: str) -> pd.DataFrame:
	try:
		df = pd.read_csv(f'{file_name}.csv').filter(regex='^(?!Unnamed|Index)')
		column_names = [col for col in df.columns]
		print(f'\nLoaded dataframe with {df.shape[0]} rows and these columns:\n\n{column_names}')
		return df
	except Exception as e:
		print(f"\nComposition list must be a csv in this script's directory:\n\n{e}")


def extract_elements(data: pd.DataFrame) -> list[str]:
	elements = [col for col in data.columns if len(col) <= 2]
	print(f'\nExample rows:\n{data.head().iloc[:, :len(elements)]}')
	return elements


def get_random_state(umap_params: dict) -> int:
	random_state = umap_params.get('random_state', None)
	if random_state is None:
		random_state = int(np.random.default_rng().integers(0, 1e9))
	return random_state

def umap_model(data: pd.DataFrame, elements: list[str], umap_params: dict):
	print('\nGenerating umap coordinates...')
	tic = time.time()
	manifold = np.array(data[elements])
	umap_model = umap.UMAP(**umap_params).fit(manifold)
	toc = time.time()
	time_taken = int(round(toc - tic, 0))
	print(f'\nLearnt Umap in {time_taken} seconds')
	return umap_model


def umap_coordinate_generator(data: pd.DataFrame, elements: list[str], umap_params: dict, umap_model) -> tuple:
	print('\nGenerating umap coordinates...')
	tic = time.time()
	manifold = np.array(data[elements])
	# umap_model = umap.UMAP(**umap_params).fit(manifold)
	coordinates = umap_model.transform(manifold)
	coordinates = pd.DataFrame(coordinates, columns=['umap_x_coordinate', 'umap_y_coordinate'])
	data_with_dra = pd.concat([data, coordinates], axis=1)
	toc = time.time()
	time_taken = int(round(toc - tic, 0))
	column_names = [col for col in data_with_dra.columns]
	print(f'\nNew dataframe with these columns in {time_taken} seconds:\n\n{column_names}')
	return data_with_dra, time_taken


def create_color_palette(elements: list) -> list[str]:
	cmap = plt.get_cmap('gist_rainbow')
	number_of_segments = len(elements)
	segment_colors = cmap(np.linspace(0, 1, number_of_segments))
	color_list = [colors.rgb2hex(segment_colors[i]) for i in range(number_of_segments)]
	return color_list


def plot_umap_data(data: pd.DataFrame, elements: list[str], color_list: list[str],
				   umap_params: dict, time_to_evaluate: int,
				   graph_title: str, file_name: str, save_image: bool):
	plt.figure(figsize=(9, 9));
	plt.axis('off')

	plt.scatter(data['umap_x_coordinate'], data['umap_y_coordinate'], color='darkgrey', edgecolor='k',
				linewidths=1.2, s=80, zorder=1)

	zorders = [i + 1 for i in range(len(elements))]
	for index, element in enumerate(elements):
		x_coordinate = data.loc[data[element] >= 0.50, 'umap_x_coordinate']
		y_coordinate = data.loc[data[element] >= 0.50, 'umap_y_coordinate']
		color = color_list[index]
		zorder = zorders[index]

		plt.scatter(x_coordinate, y_coordinate, label=element, color=color, edgecolor='k',
					linewidths=1.2, s=80, zorder=zorder)

	plt.title(f'{graph_title}', fontfamily='sans-serif',
			  fontsize=40, fontweight='heavy', ha='center', va='bottom')
	total = len(elements)
	num_columns = (total) * (total <= 6) + (int(round(total / 2, 0))) * (total > 6)
	plt.legend(loc='center', bbox_to_anchor=(0.5, -0.1), ncol=num_columns,
			   fontsize=30, handletextpad=0.2, columnspacing=0.4, markerscale=3, scatterpoints=1)
	text = [
		f'Points: {data.shape[0]}',
		f'Metric: {umap_params.get("metric")}',
		f'Seed: {umap_params.get("random_state")}',
		f'CPU Jobs: {umap_params.get("n_jobs")}',
		'',
		f'Epochs: {umap_params.get("n_epochs")}',
		f'Neighbors: {umap_params.get("n_neighbors")}',
		f'Min Dist: {umap_params.get("min_dist")}',
		f'Learn Rate: {umap_params.get("learning_rate")}',
		f'Gen Time: {time_to_evaluate}'
	]
	plt.text(1.00, 0.5, '\n'.join(text), transform=plt.gca().transAxes, fontfamily='DejaVu Sans Mono',
			 fontsize=30, ha='left', va='center')

	if save_image:
		strings = [
			f'{file_name}',
			f'points_{data.shape[0]}',
			f'jobs_{umap_params.get("n_jobs")}',
			f'epochs_{umap_params.get("n_epochs")}',
			f'nneighbors_{umap_params.get("n_neighbors")}',
			f'state_{umap_params.get("random_state")}',
			f'time_{time_to_evaluate}'
		]
		export_file_name = f'{"_".join(strings)}_96dpi2.png'
		plt.savefig('./0_UMAP_OUT.png', dpi=96, bbox_inches='tight')
		print(f'\nPlot exported as {export_file_name}.png')
	plt.show()
	plt.close()


def export_umap_data(data: pd.DataFrame, file_name: str, save_data: bool):
	if save_data:
		with open(f'{file_name}_umap.csv', 'w') as f:
			data.to_csv(f, line_terminator='\n')
			print(f'\nData exported to {file_name}_umap.csv')


def export_metadata_to_json(data: pd.DataFrame, umap_params: dict,
							time_taken: int, file_name: str, save_data: bool):
	if save_data:
		metadata = {
			'data_points': data.shape[0],
			'data_columns': list(data.columns),
			'umap_params': umap_params,
			'time_to_generate_coordinates_sec': time_taken,
			'time_to_generate_coordinates_min': int(round(time_taken / 60, 0))
		}
		with open(f"{file_name}_metadata.json", "w") as json_file:
			json.dump(metadata, json_file, indent=4)
			print(f'\nMetadata exported to {file_name}_metadata.json')


if __name__ == "__main__":
	n = 5
	file_name = f'0_{n}space'
	save_image = True
	save_data = True
	graph_title = 'UMAP'

	umap_params = {
		'n_jobs': 1,
		'random_state': 210181310,
		'n_epochs': 500,  # 400 good
		'n_components': 2,
		'n_neighbors': 150,  # 50 good
		'min_dist': 1.0,
		# 'spread': 1.0,
		'metric': 'euclidean',
		'learning_rate': 1.0
	}

	change_directory_to_current_py_file()
	data = read_data('./' + file_name)

	elements = extract_elements(data)
	color_list = create_color_palette(elements)

	random_state = get_random_state(umap_params)
	umap_params.update({'random_state': random_state})
	umap_m = umap_model(data, elements, umap_params)

	data_with_dra, time_taken = umap_coordinate_generator(data, elements, umap_params, umap_m)

	#reduced system
	mol_grid = create_mol_grid(grid_size=10, n=n)
	df = pd.DataFrame(mol_grid, columns=[str(i) for i in range(1, n + 1)])
	manifold = np.array(df[elements])
	coordinates = umap_m.transform(manifold)
	coordinates = pd.DataFrame(coordinates, columns=['umap_x_coordinate', 'umap_y_coordinate'])
	data_with_dra = pd.concat([df, coordinates], axis=1)
	plot_umap_data(data_with_dra, elements, color_list, umap_params, time_taken,
				   graph_title, file_name, save_image)
	export_umap_data(data_with_dra, file_name, save_data)
	export_metadata_to_json(df, umap_params, time_taken, file_name, save_data)
	# print(data_with_dra)
