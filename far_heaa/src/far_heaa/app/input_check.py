import streamlit as st
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))

from far_heaa.io.text_handler import TextHandler
from far_heaa.io.metadata_handler import MetadataHandler


class InputCheck:

    def __init__(self, input_str):

        self.input_str = input_str
        self.input_list = input_str.split("-")
        self.input_set = set(self.input_list)
        self.folder_path = "../database/"
        self.all_element_list = TextHandler.extract_ele_list(
            folder_path=self.folder_path, file_name="all_elements"
        )
        self.all_element_set = set(self.all_element_list)
        self.database_element_list = TextHandler.extract_ele_list(
            folder_path=self.folder_path, file_name="database_element_list"
        )
        self.database_element_set = set(self.database_element_list)

    def subset_all(self):
        return (
            self.input_set.issubset(self.all_element_set)
            and self.input_set != self.all_element_set
        )

    def subset_database(self):
        return (
            self.input_set.issubset(self.database_element_set)
            and self.input_set != self.database_element_set
        )

    def blank_element(self):
        return self.input_str[-1] == "-"

    def repetition(self):
        return len(self.input_list) != len(set(self.input_list))

    def one_element(self):
        return len(self.input_list) == 1

    def input_validation(self):

        subset_all_flag = self.subset_all()
        subset_database_flag = self.subset_database()
        blank_element_flag = self.blank_element()
        repetition_flag = self.repetition()
        one_element_flag = self.one_element()

        invalid_flag = (
            not subset_all_flag
            or not subset_database_flag
            or blank_element_flag
            or repetition_flag
            or one_element_flag
        )

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
            st.write(
                f"{', '.join(list(self.input_set.difference(self.all_element_set)))} not a valid element!"
            )
            return invalid_flag
        if not subset_database_flag:
            st.write(
                f"{', '.join(list(self.input_set.difference(self.database_element_set)))} not in database! We are working on increasing our database!"
            )
            return invalid_flag

        return invalid_flag


@st.cache_data
def get_metadata():
    mH = MetadataHandler()
    meta_data = mH.get_metadata
    return meta_data, mH
