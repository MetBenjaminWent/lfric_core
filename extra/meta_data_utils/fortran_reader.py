##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""This module contains functionality for finding and parsing LFRic diagnostic
 meta data files"""
import glob
import logging
import os
import re
from typing import Dict, Tuple

from fparser.common.readfortran import FortranFileReader
from fparser.two.Fortran2003 import Assignment_Stmt, Structure_Constructor_2, \
    Section_Subscript_List, Char_Literal_Constant, Array_Constructor, \
    Level_3_Expr, Structure_Constructor, Part_Ref
from fparser.two.parser import ParserFactory
from fparser.two.utils import walk, FparserException

from dimension_parser import parse_vertical_dimension_info, \
    translate_vertical_dimension
from entities import Field, Section, Group


class FortranMetaDataReader:
    """This encapsulates all the parsing functionality. It is give a root
    directory (Head of the repository) upon creation"""

    logger = logging.getLogger("lfric_meta_data_parser")
    vertical_dimension_definition = None
    valid_meta_data = True

    FILE_NAME_REGEX = re.compile(r"(?P<section_name>[a-z_]+?)__"
                                 r"(?P<group_name>[a-z_]+)__")

    def __init__(self, root_directory: str, meta_types_path: str):
        self.__root_dir = root_directory
        self.meta_types_path = meta_types_path
        self.meta_mod_files = None
        self.find_fortran_files()
        self.get_vertical_dimension_definition()

    def find_fortran_files(self):
        """Recursively looks for fortran files ending with "__meta_mod.f90"
        Initialises the vertical_dimension_definition variable
        Returns a list of file names"""

        self.logger.info("Scanning for fortran meta data files...")
        self.meta_mod_files = glob.glob(
            self.__root_dir + '/**/source/diagnostics_meta/**/*__meta_mod.*90',
            recursive=True)

        self.meta_mod_files.sort()
        self.logger.info("Found %i meta data files", len(self.meta_mod_files))

    def get_vertical_dimension_definition(self):
        """Looks for functions that create vertical dimension objects and reads
        the hard coded and default values"""
        self.vertical_dimension_definition = \
            parse_vertical_dimension_info(
                self.__root_dir + self.meta_types_path)

    def read_fortran_files(self) -> Tuple[Dict[str, Section], bool]:
        """Takes a list of file names (meta_mod.f90 files)
        Checks for correctness and returns the relevant fortran lines in a list
        :return Metadata: A dictionary, each key represents a fortran file and
        it's value is a list of strings, each element representing a field"""

        meta_data: Dict[str, Section] = {}
        valid_files = 0
        # Loop over each found fortran file
        for file_path in self.meta_mod_files:

            try:
                # Load the fortran file

                reader = FortranFileReader(file_path, ignore_comments=True)
                f2003_parser = ParserFactory().create(std="f2003")
                parse_tree = f2003_parser(reader)

                file_valid = True
                file_name_parts = self.FILE_NAME_REGEX.search(file_path)

                if not file_name_parts:
                    self.logger.error('Filename in path is not correct' +
                                      os.linesep + file_path)
                    self.valid_meta_data = False
                    break

                section_name = file_name_parts.group("section_name")
                group_name = file_name_parts.group("group_name")
                file_name = file_path[file_path.rfind("/")+1:]

                if section_name not in meta_data:
                    meta_data.update({section_name:
                                      Section(name=section_name)})

                group = Group(name=group_name, file_name=file_name)

                meta_data[section_name].add_group(group)

                # For every instance of a meta type object being created
                for definition in walk(parse_tree.content,
                                       types=Assignment_Stmt):

                    field, valid = self.extract_field(definition, file_name)

                    if not valid:
                        self.valid_meta_data = False
                        file_valid = False

                    if field.is_valid():
                        group.add_field(field)
                    else:
                        self.logger.error("%s is invalid. Please check",
                                          field.unique_id)
                        self.valid_meta_data = False
                        file_valid = False

                if file_valid:
                    valid_files += 1

            except FparserException as error:
                self.logger.error(": Fparser Exception %s", error)
                self.valid_meta_data = False

        if valid_files == len(self.meta_mod_files):
            self.logger.info("All %i files are valid",
                             len(self.meta_mod_files))
        else:
            self.logger.error("%i of %i files are invalid",
                              len(self.meta_mod_files) - valid_files,
                              len(self.meta_mod_files))

        return meta_data, self.valid_meta_data

    def extract_field(self, definition, file_name) -> Tuple[Field, bool]:
        """Takes an fparser object and extracts the field definition
        information
        :param definition: An fparser object that contains a field definition
        :param file_name: The name of the file that the field is declared in.
        This is needed for Field object creation
        :return Field"""

        valid_field = True
        field = Field(file_name)

        # For every instance argument in object creation
        for parameter in walk(definition, Structure_Constructor_2):

            value = None
            key = parameter.children[0].string

            try:
                # For multi line statements
                if isinstance(parameter.parent, Level_3_Expr):
                    key, value = self.extract_multi_line_statement(parameter)

                # ENUM's
                elif isinstance(parameter.children[1], Char_Literal_Constant):
                    value = parameter.children[1].string[1:-1]

                # Dimension object creation without args
                # (Structure_Constructor) or with args (Part_Ref)
                elif isinstance(parameter.children[1],
                                (Part_Ref, Structure_Constructor)):
                    value = translate_vertical_dimension(
                        self.vertical_dimension_definition,
                        parameter.children[1].string)

                # For statements with arrays in them (misc_meta_data)
                elif isinstance(parameter.children[1], Array_Constructor):
                    value = {}
                    for array in walk(parameter.children,
                                      types=Section_Subscript_List):
                        value.update({array.children[0].string[1:-1]:
                                      array.children[1].string[1:-1]})

                else:
                    # This ignore's args used in vertical dimension creation
                    key_blacklist = ["top", "bottom"]
                    if key in key_blacklist:
                        continue

                    value = parameter.children[1].string

            except Exception as error:
                if field.unique_id:
                    self.logger.warning(
                        "Key: %s on field: %s in file: %s is invalid: %s", key,
                        field.unique_id, file_name, error)
                else:
                    self.logger.warning("Key: %s in file: %s is invalid: %s ",
                                        key, file_name, error)

            # Adds the key / value to the Field object
            if hasattr(field, key):
                setattr(field, key, value)
            else:
                self.logger.error("Unexpected Field Property: %s", key)
                valid_field = False

        return field, valid_field

    @staticmethod
    def extract_multi_line_statement(statement: Level_3_Expr) -> [str, str]:
        """Accepts a fparser Level_3_Expr object as input. This object
        represents multi line statements.
        The function reassembles an arbitrary amount of lines into one
        statement
        :param statement: An fparser Level_3_Expr object
        :return key, value: The key and reassembled value as strings
        """
        value = ""
        key = None
        while isinstance(statement.parent, Level_3_Expr):
            for item in statement.parent.children:

                if isinstance(item, Structure_Constructor_2):
                    if not key:
                        key = item.children[0].string

                    # [1:-1] Removes the quotes from the string
                    value += item.children[1].string[1:-1]

                elif isinstance(item, Char_Literal_Constant):
                    value += item.children[0][1:-1]
            statement = statement.parent

        return key, value
