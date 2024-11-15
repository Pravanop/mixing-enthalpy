temp_grid = list(np.linspace(0, 2500, 40))
conv_hull = viz.grid_iterator.temp_iterator(composition=composition, temp_grid=temp_grid)

def get_mol(composition, entries):
    entries_mol = np.array([[i.composition.get_atomic_fraction(j) for j in composition] for i in entries])
    entries_mol = entries_mol[:, :-1]
    return entries_mol

temperature_dict = {}
for temperature in temp_grid:
    temp_dict = {}
    stable_entries = conv_hull[temperature].stable_entries
    stable_entries_name = [i.name for i in stable_entries]
    stable_entries_name = set(re.sub(r'\d+', '', entry) for entry in stable_entries_name)
    stable_entries_mol = get_mol(composition, stable_entries)
    stable_mol_decomp_dict = zip(stable_entries_mol, stable_entries_name)


    unique_stable_entries = list(set(stable_entries_name))

    unstable_entries = conv_hull[0].unstable_entries
    unstable_entries_decomp = [conv_hull[0].get_decomposition(i.composition) for i in unstable_entries]
    unstable_entries_mol = get_mol(composition, unstable_entries)

    decomposition_list = []
    for i in unstable_entries_decomp:
        products = list(i.keys())
        # if len(products) == 1:
        #     continue

        decomp_products = [i.name for i in products]
        decomp_products = '+'.join((re.sub(r'\d+', '', entry) for entry in decomp_products))
        decomposition_list.append(decomp_products)


    unstable_mol_decomp_dict = zip(unstable_entries_mol, decomposition_list)


    unique_decomposition_list = list(set(decomposition_list))
    total_phases = unique_decomposition_list + unique_stable_entries

    temp_dict['total_phases'] = total_phases
    temp_dict['unstable_mol_zip'] = unstable_mol_decomp_dict
    temp_dict['stable_mol_zip'] = stable_mol_decomp_dict
    temperature_dict[temperature] = temp_dict


total_unique_phases = []
for i in temperature_dict.values():
    total_unique_phases += i['total_phases']

total_unique_phases = list(set(total_unique_phases))

cmap = plt.get_cmap('viridis')  # You can choose any colormap here
colors = cmap(np.linspace(0, 1, len(total_unique_phases)))  # Generate a color for each entry

# Create a dictionary to map each entry to a color
color_dict = {entry: color for entry, color in zip(total_unique_phases, colors)}
print(color_dict)

for T, value in temperature_dict.items():
    #unstbale
    for u_mol, u_id in value['unstable_mol_zip']:
        plt.scatter(u_mol[0], T, color = color_dict[u_id])

    for s_mol, s_id in value['stable_mol_zip']:
        plt.scatter(s_mol[0], T, color= color_dict[s_id])

plt.xlim([0,1])
plt.show()