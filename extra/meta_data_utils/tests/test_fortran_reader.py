from fortran_reader import *
from diag_meta_gen import get_root_dir
import os

TEST_DATA_DIR = os.path.dirname(os.path.abspath(__file__)) + "/test_data"
meta_types_folder = "/gungho/unit-test/diagnostics_meta/meta_types/"
root_dir = get_root_dir()
test_parser = FortranMetaDataReader(root_dir, meta_types_folder)


def test_read_fortran_files_1():

    result = test_parser.read_fortran_files()
    assert "example_science_section" in result[0].keys()
    assert result[1] == True
    assert isinstance(result[0]["example_science_section"], Section)


def test_read_fortran_files_2():

    test_parser.meta_mod_files = [TEST_DATA_DIR +
                                  "/test_section__test_group__.f90"]
    result = test_parser.read_fortran_files()

    assert result[1] == True
    assert "test_section" in result[0].keys()
    assert "test_group" in result[0]["test_section"].groups
    assert "example_fields__eastward_wind" in \
           result[0]["test_section"].groups["test_group"].fields
