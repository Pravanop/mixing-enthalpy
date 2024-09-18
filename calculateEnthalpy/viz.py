from matplotlib import cm
import matplotlib as mp
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from tqdm import tqdm
import numpy as np
from calculateEnthalpy.helper_functions.grid_code import create_multinary, create_mol_grid
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram

import matplotlib.pyplot as plt

mp.rcParams['font.size'] = 12


def mol_gradation_end(member_pos, n, x):
    mols = []
    for i in x:
        subtract = np.zeros(n)
        subtract += (1 - i) / n
        subtract[member_pos] = (1 + (n - 1) * i) / n
        mols.append(subtract)
    
    return mols


def mol_gradation_end(member_pos, n, x):
    mols = []
    for i in x:
        subtract = np.zeros(n)
        subtract += (1 - i) / n
        subtract[member_pos] = (1 + (n - 1) * i) / n
        mols.append(subtract)
    
    return mols


def mol_gradation_face(member_pos, n, x):
    mol_list = []
    for mol in x:
        addition = np.zeros(n)
        addition += (n - 1 + mol) / (n * (n - 1))
        addition[member_pos] = (1 - mol) / (n)
        mol_list.append(addition)
    
    # mol_list.reverse()
    return mol_list


def assign_rolling_slice(lst, start, new_values):
    n = len(lst)
    length = len(new_values)
    
    # Loop through the new values and assign them to the original list in a rolling manner
    for i in range(length):
        lst[(start + i) % n] = new_values[i]
    
    return lst

def mol_gradation_edge(member_pos, N, n, x):
    mol_list = []
    for mol in x:
        addition = np.zeros(n)
        addition += [(1 - mol) / n] * n
        for i in member_pos:
            addition[i] = (1 / N - 1 / n) * mol + (1 / n)
        # addition = assign_rolling_slice(addition, member_pos, [(1 / N - 1 / n) * mol + (1 / n)]*N)
        mol_list.append(addition)
        
    return mol_list


def misc_temp(composition, pD, idx, x, im_list, N, flag='end'):
    if flag == 'end':
        mol_grid = mol_gradation_end(idx, len(composition), x)
    if flag == 'face':
        # mol_grid = mol_gradation_face(idx, len(composition), x)
        mol_grid = mol_gradation_edge(idx, N, len(composition), x)
    
    temp_list = []
    conv_hull = pD.make_convex_hull(composition=composition, temperature=0)
    for idx, mol in enumerate(tqdm(mol_grid)):
        misc_T = pD.find_misc_temperature(composition=composition, lattice="min", mol_ratio=mol, batch_tag=True,
                                          im=im_list,
                                          phase_flag=False, conv_hull=conv_hull)
        if isinstance(misc_T, float):
            temp_list.append(misc_T)
        else:
            temp_list.append(misc_T[0])
    
    temp_list = np.array(temp_list)
    return temp_list


def plot_line(angle_degrees, x_values, temp_list, ax, cmap, norm, angle_offset, zorder):
    angle_radians = np.radians(angle_degrees + angle_offset)
    r_values = np.array(x_values)  # Use x_values as the radius in polar coordinates
    # print(r_values[-1])
    # if flag == 'end':
    t = mp.markers.MarkerStyle(marker='|')
    t._transform = t.get_transform().rotate_deg(angle_degrees)
    # ax.scatter([angle_radians] * len(r_values), r_values, c='black', alpha=0.5, cmap=cmap, s=10, marker=t, zorder=100)
    # if flag == 'face':
    #     ax.scatter([angle_radians] * len(r_values), r_values, c=temp_list, cmap=cmap, s=140, marker = 'D', zorder = 0)
    
    for i in range(len(r_values)):
        if i + 1 <= len(r_values):
            ax.vlines(angle_radians, 0, r_values[i],
                      color=cmap(norm(temp_list[i])), linewidth=3, zorder=zorder + len(r_values) - i)


def set_composition_labels(ax, composition):
    # Set custom angular tick positions and labels
    n = len(composition)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)  # Divide the circle evenly for each composition
    off_angles = angles + np.pi / n
    off_compostions = [''.join([x for i, x in enumerate(composition) if i != idx]) for idx in range(len(composition))]
    # total = composition + off_compostions
    # total_angle = np.append(angles, off_angles)
    ax.set_xticks(off_angles)  # Set the positions of the angular ticks
    ax.set_xticklabels(off_compostions)


def label_radial_lines(ax, angles,  labels, radius=0.1):
    # Label the radial lines (inner radius)
    
    for label, angle in zip(labels, angles):
        # Use ax.text() to place labels at the specified radius and angle
        print(angle)
        ax.text(angle, radius, label, ha='center', va='center', fontsize=12, color='black')


def draw_circle_in_polar(radius, ax):
    # Create theta values from 0 to 2*pi (full circle)
    theta = np.linspace(0, 2 * np.pi, 100)
    
    # Plot the circle by setting r constant (equal to the given radius)
    ax.plot(theta, [radius] * len(theta), label=f'Circle with radius {radius}', linewidth=0.8, zorder=0, color='black',
            linestyle='-', alpha=0.1)

def find_indices(main_list, subset):
    indices = []
    for value in subset:
        try:
            index = main_list.index(value)  # Find the index of the value in the main list
            indices.append(index)
        except ValueError:
            indices.append(None)  # If the value is not found, append None or handle as needed
    return indices

def divide_circle_degrees(n):
    # Generate n equally spaced angles in degrees
    return np.linspace(0, 360, n, endpoint=False)

if __name__ == '__main__':
    
    correction = True
    equi = True
    
    if correction:
        binary_file_path = "../calcEnthalpy_old/new_phase_diagram/bokas_omegas_processed.json"
    else:
        binary_file_path = "../data/output_data/bokasCorrected_bcc_1/all_lattices_binaries.json"
    
    end_member_path = "../calcEnthalpy_old/new_phase_diagram/bokas_end_members_dict.json"
    pD = phaseDiagram(
        processed_binary_file_path=binary_file_path,
        end_member_file_path=end_member_path,
        grid_size=15,
        im_flag=False,
        correction=correction,
        equi_flag=equi)
    
    composition = ['Cr', 'V', 'Ta']
    composition.sort()
    
    n_alloy = len(composition)
    all_combs = create_multinary(element_list=composition, no_comb=list(range(2, n_alloy + 1)))
    im_list = []
    for dimensionality, alloy_list in all_combs.items():
        if pD.im_flag:
            im_list += pD.get_intermetallic(alloy_list)
    x = np.linspace(0, 1, 10)
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 10))
    cmap = cm.get_cmap('coolwarm')
    melt_T = []
    for idx, i in enumerate(composition):
        mol_ratio = [1 if idx == i else 0 for i in range(n_alloy)]
        melt_T.append(pD.tm.avg_T_melt(composition=[i], mol_ratio=mol_ratio))
    norm = Normalize(vmin=0, vmax=max(melt_T))
    # for idx, end in enumerate(composition):
    #     names = []
    #     angles = []
    n_composition = np.array(composition)
    for N in range(1, len(composition)):
        if N == 1 :
            combs = composition
        else:
            combs = list(create_multinary(element_list=composition, no_comb=[N]).values())[0]
        
        angles = divide_circle_degrees(len(combs))
        for idx, i in enumerate(combs):
            # angle = idx*360/n_alloy + (N-1)*360/len(combs)
            angle = angles[idx] + (N-1)*10
            print(angle, i)
            # print(i)
            temp_i = i.split('-')
            # print(temp_i)
            member_pos = find_indices(composition, temp_i)
            # print(member_pos)
            angle_offset = 0
            misc_temp_list = misc_temp(composition, pD, member_pos, x, im_list, flag='face', N=N)
            plot_line(angle_degrees=angle, x_values=x*N*1.5, temp_list=misc_temp_list, ax=ax, cmap=cmap, norm=norm,
                      angle_offset=angle_offset, zorder=1)
            angle_radians = np.radians(angle + angle_offset)
            if angle > 90 and angle < 270:
                rotation = angle + 180  # Flip the text for the left side of the circle
            else:
                rotation = angle
            ax.text(angle_radians, N*1.5*(1+0.2), i, ha='center', va='center', color='black', rotation=rotation)
            
            # angles.append(angle)
            
        # print(angles)
        # label_radial_lines(ax, angles = angles, labels = combs, radius=N)
    
    for N in range(1, len(composition)):
        draw_circle_in_polar(radius=N*1.5, ax=ax)
        
    
    ax.scatter(0, 0, c=cmap(norm(misc_temp_list[0])), cmap=cmap, marker='o', s=140, zorder=100, edgecolor='black')
    
    # draw_circle_in_polar(radius=2.0, ax=ax)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines['polar'].set_visible(False)
    # label_radial_lines(ax, composition, radius=1.2)
    # set_composition_labels(ax, composition)
    ax.grid(False)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # We need this for colorbar to work
    cbar = plt.colorbar(sm, ax=ax, aspect=30, fraction=0.05, orientation='horizontal')  # Unified colorbar
    cbar.set_label('$T_{misc}$', fontsize=12)
    plt.savefig(f'./{"".join(composition)}.png', dpi=300)
