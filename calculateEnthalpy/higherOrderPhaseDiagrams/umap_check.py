import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import scipy
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter
from tqdm import tqdm
import mpltern
from calculateEnthalpy.helper_functions.grid_code import create_multinary, create_mol_grid
from calculateEnthalpy.helper_functions.phase_diagram import phaseDiagram




if __name__ == '__main__':
    composition = ['Fe', 'Ta', 'W']

    n_alloy = len(composition)
    all_combs = create_multinary(element_list=composition, no_comb=list(range(2, n_alloy + 1)))
    im_list = []

    mol_grid = create_mol_grid(3, 15)
    mol_grid = np.round(np.array(mol_grid),4)

    correction = True
    equi = False

    if correction:
        binary_file_path = "../../calcEnthalpy_old/new_phase_diagram/bokas_omegas_processed.json"
    else:
        binary_file_path = "../../data/output_data/bokasCorrected_bcc_1/all_lattices_binaries.json"

    end_member_path = "../../calcEnthalpy_old/new_phase_diagram/bokas_end_members_dict.json"
    pD = phaseDiagram(
        processed_binary_file_path=binary_file_path,
        end_member_file_path=end_member_path,
        grid_size=10,
        im_flag=True,
        correction=correction,
        equi_flag=equi)
    for dimensionality, alloy_list in all_combs.items():
        if pD.im_flag:
            im_list += pD.get_intermetallic(alloy_list)

    temp_list = []
    for mol in tqdm(mol_grid, desc='Iterating',smoothing=0.8):
        misc_T = pD.find_misc_temperature(composition=composition, lattice="min", mol_ratio=mol, batch_tag=True, im=im_list, phase_flag=False)
        if isinstance(misc_T, float):
            temp_list.append(misc_T)
        else:
            temp_list.append(misc_T[0])

    tm = temp_list

    t, l, r = mol_grid[:, 0], mol_grid[:, 1], mol_grid[:, 2]
    t = t.reshape(-1, 1)
    l = l.reshape(-1, 1)
    r = r.reshape(-1, 1)
    tm = np.array(tm).reshape(-1, 1)
    fig = plt.figure()

    ax = fig.add_subplot(projection="ternary")
    ax.grid()
    cax = ax.inset_axes([1.03, 0.1, 0.05, 0.9], transform=ax.transAxes)
    pc = ax.scatter(t, l, r, c = tm, cmap = "coolwarm")
    data = np.concatenate([t, l, r, tm], axis = 1)
    # # data = data.reshape((int(data.size/4),4))
    df = pd.DataFrame(data, columns = ["t", "l", "r", "tm"])
    df.to_csv("extract.csv")
    # # data = scipy.ndimage.zoom(data, 2)
    # sigma=0.15
    # data = gaussian_filter(data, sigma)
    # ps = ax.tricontour(t,l,r,tm, cmap = "coolwarm")
    # ps = ax.tricontour(data[:,0], data[:,1], data[:, 2], data[:,3], cmap = "coolwarm")
    colorbar = fig.colorbar(pc, cax=cax)
    colorbar.set_label('T_misc', rotation=270, va='baseline')
    ax.set_tlabel(f"{composition[0]}")
    ax.set_llabel(f"{composition[1]}")
    ax.set_rlabel(f"{composition[2]}")

    plt.savefig(f"../../plots/phase_diagrams/{'-'.join(composition)}")
# plt.show()
