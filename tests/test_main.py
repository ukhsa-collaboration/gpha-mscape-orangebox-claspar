import argparse
from unittest.mock import patch

from claspar import main

PATH_TO_TEST_TABLE = "tests/test_profile_tables/test_profile_tables.xlsx"


@patch(
    "argparse.ArgumentParser.parse_args",
    return_value=argparse.Namespace(
        sample_id="ID-12345678",
        output_dir="tests/end_to_end",
        server="mscape",
        samplesheet_path="tests/test_samplesheet_main_input.tsv",
        profile_table_spreadsheet_path=PATH_TO_TEST_TABLE,
        log_file=None,
        config=None,
        database_path=None,
    ),
)
def test_end_to_end(mock_args, capsys):
    exitcode = main.main()
    assert exitcode == 0
