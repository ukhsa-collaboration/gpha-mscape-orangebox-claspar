import argparse
import json
from pathlib import Path
from unittest.mock import patch

from claspar import main

root = Path(__file__).parents[1]

PATH_TO_TEST_TABLE = Path(root / "tests/test_profile_tables/profile_202303.0.0.json")
PATH_TO_TEST_DATA = Path(root / "tests/test_unit_test_metadata.json")

with PATH_TO_TEST_DATA.open("r") as metadata:
    MOCK_ONYX_RECORD = json.load(metadata)


def test_taxaplease():
    from taxaplease import TaxaPlease

    tp = TaxaPlease()
    print("Your taxaplease database version is: ")
    print(tp.get_current_taxonomy_url_from_database())


@patch("onyx_analysis_helper.onyx_analysis_helper_functions.OnyxClient.get")
@patch(
    "argparse.ArgumentParser.parse_args",
    return_value=argparse.Namespace(
        sample_id="ID-12345678",
        output_dir="tests/end_to_end",
        server="mscape",
        samplesheet_path=None,
        profile_json_path=PATH_TO_TEST_TABLE,  # patch in the variable not the arg
        log_file=None,
        config=None,
        database_path=None,
    ),
)
def test_end_to_end(mock_args, mock_query, capsys, caplog):
    # First mock the onyx query results.
    mock_query.return_value = MOCK_ONYX_RECORD

    exitcode = main.main()
    print(caplog.text)
    assert exitcode == 0
