from json_meta_data import *

TEST_DATA_DIR = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"


def test_write_json_meta():
    test_data = {'section_name': 'test_section', 'unique_id': 'test_id'}

    write_json_meta(test_data, TEST_DATA_DIR, )
    with open(os.path.join(TEST_DATA_DIR, "rose.meta")) as json_test_data:
        test_data = json.load(json_test_data)
        assert test_data == {
            'meta_data': {
                'section_name': 'test_section',
                'unique_id': 'test_id'
            },
            'checksum': 'md5: 1461c4f998be95df7f0512b9ed772ce0'}


def test_set_checksum():
    test_dict = {}
    test_dict_2 = {"a key": "a value"}

    set_checksum(test_dict)
    set_checksum(test_dict_2)

    assert test_dict['checksum'] == "md5: 99914b932bd37a50b983c5e7c90ae93b"
    assert test_dict_2['checksum'] == "md5: 9eae1793e5a5fbd89e0b737c67173a46"
