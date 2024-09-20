import os

import numpy as np
from matplotlib import pyplot as plt, cm
from matplotlib.colors import Normalize

from calcEnthalpy_package.grids_and_combinations.combination_generation import MultinaryCombinations
from calcEnthalpy_package.grids_and_combinations.grid_creation import CompositionGrid
from calcEnthalpy_package.io.json_handler import JSONHandler
from calcEnthalpy_package.math_operations.thermo_calculations import ThermoMaths
from calcEnthalpy_package.phase_diagram.grid_iterators import GridIterator
from calcEnthalpy_package.math_operations.polar_calculations import PolarMaths as pm


class PolarVisualizations:
    
    def __init__(self,
                 composition,
                 lattice,
                 meta_data,
                 save_flag,
                 ):
        
        self.composition = composition
        self.lattice = lattice
        self.meta_data = meta_data
        self.save_flag = save_flag
        
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
        
        self.cmap = cm.get_cmap('coolwarm')
        
        t_max = max([self.tm.avg_T_melt(i, mol_ratio=[]) for i in self.composition])
        self.temp_grid = list(np.linspace(0, t_max, 30))
        self.norm = Normalize(vmin=0, vmax=t_max)
        self.conv_hull = self.grid_iterator.temp_iterator(composition=self.composition,
                                                          temp_grid=self.temp_grid)
    
    @staticmethod
    def assign_rolling_slice(lst, start, new_values):
        n = len(lst)
        length = len(new_values)
        
        # Loop through the new values and assign them to the original list in a rolling manner
        for i in range(length):
            lst[(start + i) % n] = new_values[i]
        
        return lst
    
    def misc_temp(self,
                  member_pos,
                  x,
                  N,
                  flag):
        n = len(self.composition)
        
        if flag == 'add':
            mol_grid = CompositionGrid.create_high_sym_mol_grid(x=x,
                                                                n=n,
                                                                N=N,
                                                                change_idx=member_pos)
        if flag == 'transmutate':
            mol_grid = CompositionGrid.create_mol_grid_transmutation(x=x,
                                                                     n=n,
                                                                     transmutation_indice=member_pos)
        
        mol_grid, misc_temp = self.grid_iterator.misc_temperature_across_grid(
            composition=self.composition,
            mol_grid_size=mol_grid,
            lattice=self.lattice,
            phase_flag=False,
            conv_hull=self.conv_hull,
            temp_grid=self.temp_grid,
        )
        return misc_temp
    
    def plot_line(self, angle_degrees, x_values, temp_list, ax, cmap, norm, zorder):
        angle_radians = np.radians(angle_degrees)
        r_values = np.array(x_values)
        for i in range(len(r_values)):
            if i + 1 <= len(r_values):
                ax.vlines(angle_radians,
                          0,
                          r_values[i],
                          color=cmap(norm(temp_list[i])),
                          linewidth=7,
                          zorder=zorder + len(r_values) - i)
    
    @staticmethod
    def distance_calculator(n, N):
        return 2 * 1 / np.sqrt((n - N))
    
    @staticmethod
    def draw_circle_in_polar(radius, ax):
        # Create theta values from 0 to 2*pi (full circle)
        theta = np.linspace(0, 2 * np.pi, 100)
        
        # Plot the circle by setting r constant (equal to the given radius)
        ax.plot(theta, [radius] * len(theta), label=f'Circle with radius {radius}', linewidth=0.8, zorder=0,
                color='black',
                linestyle='--', alpha=0.1)
    
    @staticmethod
    def find_indices(main_list, subset):
        indices = []
        for value in subset:
            try:
                index = main_list.index(value)  # Find the index of the value in the main list
                indices.append(index)
            except ValueError:
                indices.append(None)  # If the value is not found, append None or handle as needed
        return indices
    
    @staticmethod
    def get_n_colors_from_cmap(cmap_name, N):
        cmap = plt.get_cmap(cmap_name)  # Get the colormap
        colors = [cmap(i) for i in np.linspace(0, 1, N)]  # Get N evenly spaced colors
        return colors
    
    @staticmethod
    def total_num(n):
        return 2 ** n - 2
    
    def total_plot(self):
        n_alloy = len(self.composition)
        x = np.linspace(0, 1, 20)
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 10))
        melt_T = []
        for idx, i in enumerate(self.composition):
            mol_ratio = [1 if idx == i else 0 for i in range(n_alloy)]
            melt_T.append(self.tm.avg_T_melt(composition=[i], mol_ratio=mol_ratio))
        norm = Normalize(vmin=0, vmax=max(melt_T))
        
        line_colors = self.get_n_colors_from_cmap('Dark2', n_alloy - 1)
        angles = pm.divide_circle_degrees(self.total_num(n_alloy))
        count = 0
        for idx2, N in enumerate(range(1, len(self.composition))):
            if N == 1:
                combs = self.composition
            else:
                combs = list(MultinaryCombinations.create_multinary(element_list=self.composition,
                                                                    no_comb=[N]).values())[0]
            
            for idx, i in enumerate(combs):
                angle = angles[count]
                temp_i = i.split('-')
                member_pos = self.find_indices(self.composition, temp_i)
                angle_offset = 0
                misc_temp_list = self.misc_temp(member_pos=member_pos,
                                                x=x,
                                                flag='add',
                                                N=N)
                self.plot_line(angle_degrees=angle, x_values=x * self.distance_calculator(n_alloy, N),
                               temp_list=misc_temp_list,
                               ax=ax,
                               cmap=self.cmap,
                               norm=norm,
                               zorder=1)
                angle_radians = np.radians(angle + angle_offset)
                
                ax.vlines(angle_radians, 0, (x * self.distance_calculator(n_alloy, n_alloy - 1))[-1], linestyles='-',
                          color=line_colors[idx2], zorder=0, alpha=0.7, linewidth=1.2)
                if 90 < angle < 270:
                    rotation = angle + 180  # Flip the text for the left side of the circle
                else:
                    rotation = angle
                
                ax.text(angle_radians, self.distance_calculator(n_alloy, n_alloy - 1) + N * 0.1, i, ha='center',
                        va='center',
                        color='black', rotation=rotation)
                
                count += 1
        
        for N in range(1, len(self.composition)):
            self.draw_circle_in_polar(radius=self.distance_calculator(n_alloy, N), ax=ax)
        
        ax.scatter(0, 0, color=self.cmap(norm(misc_temp_list[0])), marker='o', s=140, zorder=100, edgecolor='black')
        
        ax.set_yticks([])
        ax.set_xticks([])
        ax.spines['polar'].set_visible(False)
        ax.grid(False)
        sm = plt.cm.ScalarMappable(cmap=self.cmap, norm=norm)
        sm.set_array([])  # We need this for colorbar to work
        cbar = plt.colorbar(sm, ax=ax, aspect=30, fraction=0.05, orientation='horizontal')  # Unified colorbar
        cbar.set_label('$T_{misc}$', fontsize=12)
        if self.save_flag:
            if not os.path.exists("../plots/polar_plots"):
                os.mkdir("../plots/polar_plots")
            plt.savefig(f'../plots/polar_plots/{"".join(sorted(self.composition))}_total.png', dpi=300)
        
        return ax, fig
    
    def plot_colored_secant(self, ax, radius, angle_start, angle_end, misc_temp_list_sec, cmap, norm, zorder=50):
        n_segments = len(misc_temp_list_sec)
        
        x_start, y_start = pm.polar_to_cartesian(radius, angle_start)
        x_end, y_end = pm.polar_to_cartesian(radius, angle_end)
        total_x = np.linspace(x_start, x_end, n_segments)
        slope = (y_end - y_start) / (x_end - x_start)
        y = slope * (total_x - x_start) + y_start
        radius, theta = pm.cartesian_to_polar(total_x, y)
        
        for idx in range(len(misc_temp_list_sec)):
            
            if idx + 1 < len(misc_temp_list_sec):
                x_coord = [theta[idx], theta[idx + 1]]
                y_coord = [radius[idx], radius[idx + 1]]
                ax.plot(x_coord, y_coord, c=cmap(norm(misc_temp_list_sec[idx])), linewidth=7,
                        zorder=zorder + len(misc_temp_list_sec) - idx)
    
    def plot_N(self, N_ind, transmute_indices):
        n_alloy = len(self.composition)
        x = np.linspace(0, 1, 20)
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 10))
        
        melt_T = []
        for idx, i in enumerate(self.composition):
            mol_ratio = [1 if idx == i else 0 for i in range(n_alloy)]
            melt_T.append(self.tm.avg_T_melt(composition=[i], mol_ratio=mol_ratio))
        
        line_colors = self.get_n_colors_from_cmap('Dark2', n_alloy - 1)
        
        count = 0
        if N_ind == 1:
            combs = self.composition
        else:
            combs = list(MultinaryCombinations.create_multinary(element_list=self.composition,
                                                                no_comb=[N_ind],
                                                                sort=False).values())[0]
        combs_n_N = []
        for idx, i in enumerate(combs):
            i = i.split('-')
            leftover = list(set(self.composition).difference(set(i)))
            
            combs_n_N.append('-'.join(leftover))
        
        for idx2, comb in enumerate([combs, combs_n_N]):
            
            if len(comb) % 2 == 0:
                angles = pm.divide_circle_degrees(len(comb)) / 2
            else:
                angles = pm.divide_circle_degrees(len(comb))
            
            for idx, i in enumerate(comb):
                angle = angles[idx] + idx2 * 180
                temp_i = i.split('-')
                member_pos = self.find_indices(self.composition, temp_i)
                N = len(temp_i)
                misc_temp_list = self.misc_temp(member_pos=member_pos,
                                                x=x,
                                                flag='add',
                                                N=N)
                self.plot_line(angle_degrees=angle,
                               x_values=x * self.distance_calculator(n_alloy, N),
                               temp_list=misc_temp_list,
                               ax=ax,
                               cmap=self.cmap,
                               norm=self.norm,
                               zorder=1)
                angle_radians = np.radians(angle)
                
                ax.vlines(angle_radians, 0, (x * self.distance_calculator(n_alloy, n_alloy - 1))[-1], linestyles='-',
                          color=line_colors[idx2], zorder=0, alpha=0.7, linewidth=1.2)
                if angle > 90 and angle < 270:
                    rotation = angle + 180  # Flip the text for the left side of the circle
                else:
                    rotation = angle
                
                ax.text(angle_radians, self.distance_calculator(n_alloy, n_alloy - 1) + N * 0.1, i, ha='center',
                        va='center',
                        color='black', rotation=rotation)
                
                #
                if transmute_indices:
                    if idx2 != 0 and idx == min(transmute_indices):
                        # Calculate the secondary miscibility temperature list for plotting
                        misc_temp_list_sec = self.misc_temp(member_pos=transmute_indices,
                                                            x=x,
                                                            flag='transmutate',
                                                            N=N)
                        
                        angle_start = angle_radians
                        angle_end = np.radians(angles[transmute_indices[1]] + idx2 * 180)
                        
                        radius = (x * self.distance_calculator(n_alloy, n_alloy - 1))[-1]
                        
                        # Plot the colored secant line
                        self.plot_colored_secant(ax, radius, angle_start, angle_end, misc_temp_list_sec, cmap=self.cmap,
                                                 norm=self.norm)
                #
                
                count += 1
        
        ax.scatter(0, 0, color=self.cmap(self.norm(misc_temp_list[0])), marker='o', s=140, zorder=100, edgecolor='black')
        
        ax.set_yticks([])
        ax.set_xticks([])
        ax.spines['polar'].set_visible(False)
        ax.grid(False)
        sm = plt.cm.ScalarMappable(cmap=self.cmap, norm=self.norm)
        sm.set_array([])  # We need this for colorbar to work
        cbar = plt.colorbar(sm, ax=ax, aspect=30, fraction=0.05, orientation='horizontal')  # Unified colorbar
        cbar.set_label('$T_{misc}$', fontsize=12)
        if self.save_flag:
            if not os.path.exists("../plots/polar_plots"):
                os.mkdir("../plots/polar_plots")
            plt.savefig(f'../plots/polar_plots/{"".join(sorted(self.composition))}_{N_ind}.png', dpi=300)
        
        return ax, fig
