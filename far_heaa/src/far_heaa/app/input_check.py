import streamlit as st
import sys
import os
sys.path.insert(0, os.path.abspath('../..'))


from far_heaa.io.text_handler import TextHandler
from far_heaa.io.metadata_handler import MetadataHandler


def input_validation(input_str):
	folder_path = '../database/'
	input_list = input_str.split('-')
	input_set = set(input_list)
	
	all_element_list = TextHandler.extract_ele_list(folder_path=folder_path, file_name='all_elements')
	all_element_set = set(all_element_list)
	
	database_element_list = TextHandler.extract_ele_list(folder_path=folder_path,
														 file_name='database_element_list')
	database_element_set = set(database_element_list)
	
	subset_all_flag = input_set.issubset(all_element_set) and input_set != all_element_set
	subset_database_flag = input_set.issubset(database_element_set) and input_set != database_element_set
	blank_element_flag = input_str[-1] == '-'
	repetition_flag = len(input_list) != len(set(input_list))
	one_element_flag = len(input_list) == 1
	invalid_flag = not subset_all_flag or not subset_database_flag or blank_element_flag or repetition_flag or one_element_flag
	
	if one_element_flag:
		st.write("Please provide more elements!")
		return invalid_flag
	if blank_element_flag:
		st.write("Blank is not an element!")
		return invalid_flag
	if repetition_flag:
		st.write("Repeated element!")
		return invalid_flag
	if not subset_all_flag:
		st.write(f"{', '.join(list(input_set.difference(all_element_set)))} not a valid element!")
		return invalid_flag
	if not subset_database_flag:
		st.write(
			f"{', '.join(list(input_set.difference(database_element_set)))} not in database! We are working on increasing our database!")
		return invalid_flag
	
	return invalid_flag, input_list


@st.cache_data
def get_metadata():
	mH = MetadataHandler()
	meta_data = mH.get_metadata
	return meta_data
