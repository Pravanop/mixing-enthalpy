from calculateEnthalpy.helper_functions.data_utils import DataUtils


def enthalpy_model_lookup(source, lattice, folder_path):
    binary_dict = DataUtils.load_json(folder_path=folder_path, lattice=lattice, source=source)
    extract_ele_lists = DataUtils.extract_ele_list(folder_path=folder_path, lattice=lattice, source=source)
    ele_dict = {}
    for idx, i in enumerate(extract_ele_lists):
        binary_dict.update({
            f'{i}-{i}': 0
        })
        ele_dict[i] = idx + 1
    return binary_dict, ele_dict




