
from dimension_parser import *
from diag_meta_gen import get_root_dir
import os

TEST_DIR = os.path.dirname(os.path.abspath(__file__)) + "/test_data"
enum_test_file = "/enum_test_file"
dimension_test_file = "/dimension_parser_test_file"


def test_read_enum():
    assert read_enum(TEST_DIR + enum_test_file) == \
           ["ONE", "TWO"]


def test_get_default_values():
    test_file = "/get_default_values_test_file"
    path = get_root_dir() + "/gungho/unit-test/diagnostics_meta/meta_types/"
    levels = read_enum(path + "levels_enum_mod.f90")
    reader = FortranFileReader(TEST_DIR + test_file)
    parse_tree = F2003_PARSER(reader)

    results = []

    for statement in walk(parse_tree.content, If_Construct):
        results.append(get_default_values(statement, path, levels))

    assert results == [('test_value_1', 'SOME_TEST_VALUE'),
                       ('Test_value_2', 'SOME_OTHER_TEST_VALUE')]


def test_get_hard_coded_values():
    test_file = "/get_hard_coded_value_test_file"
    reader = FortranFileReader(TEST_DIR + test_file)

    parse_tree = F2003_PARSER(reader)

    results = []

    for statement in walk(parse_tree.content, Assignment_Stmt):
        results.append(get_hard_coded_values(statement))
    assert results == [{'another_test_string': '',
                        'test_name': 'test_string',
                        'test_symbol': 'TEST_SYMBOL'}]


def test_parse_vertical_dimension_info():
    root_dir = get_root_dir() + "/gungho/unit-test/diagnostics_meta/meta_types/"
    dimension_def = parse_vertical_dimension_info(root_dir)

    assert "model_depth_dimension" in dimension_def.keys()
    assert "model_height_dimension" in dimension_def.keys()
    assert "fixed_depth_dimension" in dimension_def.keys()
    assert "fixed_height_dimension" in dimension_def.keys()


def test_translate_vertical_dimension():
    root_dir = get_root_dir() + "/gungho/unit-test/diagnostics_meta/meta_types/"
    dimension_def = parse_vertical_dimension_info(root_dir)

    model_height = translate_vertical_dimension(dimension_def,
                                                "model_height_dimension()")
    model_depth = translate_vertical_dimension(dimension_def,
                                                "model_depth_dimension()")
    fixed_height = translate_vertical_dimension(dimension_def,
                                                "fixed_height_dimension()")
    fixed_depth = translate_vertical_dimension(dimension_def,
                                                "fixed_depth_dimension()")

    assert model_height == {"standard_name": "height",
                            "units": "m",
                            "top_arg": "TOP_ATMOSPHERIC_LEVEL",
                            "bottom_arg": "BOTTOM_ATMOSPHERIC_LEVEL",
                            "positive": "POSITIVE_UP"}
    assert model_depth == {"standard_name": "depth",
                            "units": "m",
                            "top_arg": "TOP_SOIL_LEVEL",
                            "bottom_arg": "BOTTOM_SOIL_LEVEL",
                            "positive": "POSITIVE_DOWN"}
    assert fixed_height == {"standard_name": "height",
                            "units": "m",
                            "positive": "POSITIVE_UP"}
    assert fixed_depth == {"standard_name": "depth",
                            "units": "m",
                            "positive": "POSITIVE_DOWN"}
