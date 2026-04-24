import argparse
from unittest.mock import patch

from claspar import main

PATH_TO_TEST_TABLE = "tests/test_profile_tables/test_profile_tables.xlsx"

MOCK_ONYX_RECORD: dict[str, str] = {
    "climb-id": "ID-123456",
    "site": "test",
    "published_date": "2026-01-01",
    "classifier_version": "1.0.0",
    "classifier_db_date": "1970-01-01",
    "ncbi_taxonomy_date": "1970-01-01",
    "scylla_version": "1.0.0",
    "sylph_db_version": "1.0.0",
    "alignment_db_version": "1.0.0",
}


@patch("onyx_analysis_helper.onyx_analysis_helper_functions.OnyxClient.get")
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
def test_end_to_end_with_samplesheet(mock_args, mock_query, capsys, caplog):
    # First mock the onyx query results.
    mock_query.return_value = MOCK_ONYX_RECORD

    exitcode = main.main()
    print(caplog.text)
    assert exitcode == 0
