import os
import json
import yaml
import shutil

def load_yaml_to_dict(yaml_file_path):
    try:
        with open(yaml_file_path, 'r') as file:
            data = yaml.safe_load(file)
        return data
    except FileNotFoundError:
        print(f"Error: The file {yaml_file_path} was not found.")
    except yaml.YAMLError as e:
        print(f"Error parsing YAML file: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

def load_json_to_dict(json_file_path):
    try:
        with open(json_file_path, 'r') as file:
            data = json.load(file)
        return data
    except json.JSONDecodeError:
        print("Error: The file is not a valid JSON.")
    except FileNotFoundError:
        print(f"Error: The file {json_file_path} was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def concatenate_files(file_paths, output_file_path):
    try:
        with open(output_file_path, 'w') as output_file:
            for file_path in file_paths:
                with open(file_path, 'r') as file:
                    output_file.write(file.read())
    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def get_element_symbols(structure):
    element_symbols = [site.specie.symbol for site in structure]
    return list(set(element_symbols))

def get_filename_without_extension(file_path) :
	base_name = os.path.basename(file_path)
	file_name , _ = os.path.splitext(base_name)
	return file_name

def create_directory(directory_path) :
	if not os.path.exists(directory_path) :
		os.makedirs(directory_path)
	else :
		print(f"Directory already exists: {directory_path}")

def string_to_file(text, file_path):
    with open(file_path, 'w') as file:
        file.write(text)

def copy_file(src_file, dst_file):
    shutil.copy2(src_file , dst_file)
    
def identify_file(path):
    if os.path.isdir(path):
        return False
    elif os.path.isfile(path):
        return True
    else:
        return "Unknown or does not exist"