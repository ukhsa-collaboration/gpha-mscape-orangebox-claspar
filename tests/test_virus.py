import json

import pandas as pd
import pytest
from taxaplease import TaxaPlease

from claspar import virus

tp = TaxaPlease()

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)


@pytest.fixture(scope="module")
def metadata_json():
    with open("tests/test_unit_test_metadata.json") as json_file:  # noqa: PTH123
        metadata = json.loads(json_file.read())
    assert metadata["alignment_results"]
    return metadata


def test_input(metadata_json):
    print(metadata_json["alignment_results"])


class TestViralParser:
    @pytest.fixture(autouse=True)
    def viral_aligner_test_data(self, metadata_json):
        self.test_input_data = pd.DataFrame(metadata_json["alignment_results"])

    @pytest.fixture(autouse=True)
    def viral_aligner_test_thresholds(self):
        thresholds = {
            "EVENNESS_VALUE": 25,
            "COVERAGE_1X": 25,
            "UNIQUELY_MAPPED_READS": 0,
            "MEAN_READ_IDENTITY": 80,
            "MEAN_ALIGNMENT_LENGTH": 500,
        }
        self.thresholds = thresholds

    @pytest.fixture(autouse=True)
    def test_instance_1(self, viral_aligner_test_thresholds, viral_aligner_test_data):
        self.test_class_instance = virus.VirusClasPar(
            sample_id="ID-12345678",
            original_viral_aligner_df=self.test_input_data,
            virus_thresholds_dict=self.thresholds,  # type: ignore
            profiles_dict={
                "ProfileX": {
                    3428501: {"taxon": "Enterovirus alpharhino", "rank": "species"},
                    3241406: {"taxon": "Mastadenovirus blackbeardi", "rank": "species"},
                },
                "ProfileY": {
                    10509: {"taxon": "Mastadenovirus", "rank": "genus"},
                },
                "ProfileZ": {
                    2731619: {"taxon": "Caudoviricetes", "rank": "class"},
                },
            },
            taxaplease_instance=tp,
            server="mscape",
        )

    def test_test_class(self):
        assert (s := self.test_class_instance.sample_id) == "ID-12345678", (
            f"Sanity checking the test class setup failed - expected ID-12345678, got {s}"
        )

    def test_filter_viral_aligner(self):
        va_processed_df = self.test_class_instance._process_viral_aligner()
        assert (s := self.test_class_instance.viral_aligner_df.shape[0]) == 10, (
            f"Expected all results to pass the filter, got {s}"
        )
        assert (a := va_processed_df.shape) == (10, 18), f"Expected 35 rows, got {a}"
        assert all(col in va_processed_df.columns.values for col in ["taxon_rank", "confidence"])

    def test_get_viral_aligner_results(self):
        headline, results, va_processed_df = self.test_class_instance._get_viral_aligner_results()
        print(headline, results)
        expected_headline = (
            "Viral Aligner classified 10 viral taxa; ProfileX-high, ProfileY-high, ProfileZ-high, ProfileZ-low, "
            "Unknown-high"
        )
        assert expected_headline in headline

        assert all(len(inner.keys()) == 4 for inner in results.values())
        assert len(results) == 4

    def test_get_virus_analysis_table(self):
        analysis_table = self.test_class_instance.get_virus_analysis_table()
        assert (p := analysis_table.pipeline_name) == "ClasPar", f'Expected pipeline name "ClasPar", got "{p}"'
        assert (n := analysis_table.name) == "claspar-viralaligner-virus", (
            f'Expected name "virus-classifier-parser", got {n}'
        )

    def test_save_outputs_to_csv(self, tmp_path):
        self.test_class_instance.save_outputs_to_csv(tmp_path)
        print(f"Saving to {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.test_class_instance.get_virus_analysis_table()
        filename = tmp_path / f"{self.test_class_instance.sample_id}_viral_aligner_analysis_fields.json"
        self.test_class_instance.analysis_table.write_analysis_to_json(filename)
        print(f"Saving json to {filename}")


class TestNoVirus:
    @pytest.fixture(autouse=True)
    def viral_aligner_test_thresholds(self):
        thresholds = {
            "EVENNESS_VALUE": 25,
            "COVERAGE_1X": 25,
            "UNIQUELY_MAPPED_READS": 0,
            "MEAN_READ_IDENTITY": 80,
            "MEAN_ALIGNMENT_LENGTH": 500,
        }
        self.thresholds = thresholds

    @pytest.fixture(autouse=True)
    def test_instance_1(self, viral_aligner_test_thresholds):
        self.test_no_data_instance = virus.VirusClasPar(
            sample_id="ID-87654321",
            original_viral_aligner_df=pd.DataFrame(),
            virus_thresholds_dict=self.thresholds,  # ty:ignore[invalid-argument-type]
            profiles_dict={},
            taxaplease_instance=tp,
            server="mscape",
        )

    def test_get_results_with_no_results(self):
        # If the alignment has been run but there are no results, or it has not been run, both return an empty dataframe
        headline, results, va_processed_df = self.test_no_data_instance._get_viral_aligner_results()
        assert "Viral Aligner classified 0 viral taxa." in headline
        assert results == {}, f"Expected empty dict, got {results}"
        print(f"\nHeadline is '{headline}'\n and results are: {results}")

    def test_filter_viral_aligner_no_results(self):
        filtered_df = self.test_no_data_instance._process_viral_aligner()
        assert filtered_df.empty, f"Expected empty dataframe, got {filtered_df}"
        print(f"\nFunction filter_viral_aligner with empty dataframe correctly gives: {filtered_df}.")

    def test_save_outputs_to_csv(self, tmp_path):
        self.test_no_data_instance.save_outputs_to_csv(tmp_path)
        print(f"Saving broken data csvs {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.test_no_data_instance.get_virus_analysis_table()

        filename = tmp_path / f"{self.test_no_data_instance.sample_id}_no_viral_aligner_analysis_fields.json"
        self.test_no_data_instance.analysis_table.write_analysis_to_json(filename)
        print(f"Saving json to {filename}")


class TestBrokenInput:
    @pytest.fixture(autouse=True)
    def viral_aligner_test_thresholds(self):
        thresholds = {
            "EVENNESS_VALUE": 25,
            "COVERAGE_1X": 25,
            "UNIQUELY_MAPPED_READS": 0,
            "MEAN_READ_IDENTITY": 90,
            "MEAN_ALIGNMENT_LENGTH": 500,
        }
        self.thresholds = thresholds

    @pytest.fixture(autouse=True)
    def viral_aligner_broken_test_data(self, metadata_json):
        viral_aligner_test_data = pd.DataFrame(metadata_json["alignment_results"])
        broken_test_data = viral_aligner_test_data.drop("mean_alignment_length", axis=1)
        self.broken_input_data = broken_test_data

    @pytest.fixture(autouse=True)
    def test_instance_1(self, viral_aligner_test_thresholds, viral_aligner_broken_test_data):
        self.test_broken_data_instance = virus.VirusClasPar(
            sample_id="ID-9999999",
            original_viral_aligner_df=self.broken_input_data,
            virus_thresholds_dict=self.thresholds,  # type: ignore
            profiles_dict={},
            taxaplease_instance=tp,
            server="mscape",
        )

    def test_test_class(self):
        assert (s := self.test_broken_data_instance.sample_id) == "ID-9999999", (
            f"Sanity checking the test class setup failed - expected ID-9999999, got {s}"
        )

    def test_nothing_breaks_if_missing_column(self):
        assert (r := self.test_broken_data_instance.results) == {}, f"Expected empty results dict, got {r}"
        self.test_broken_data_instance.get_virus_analysis_table()

    def test_save_outputs_to_csv(self, tmp_path):
        self.test_broken_data_instance.save_outputs_to_csv(tmp_path)
        print(f"Saving broken data csvs {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.test_broken_data_instance.get_virus_analysis_table()

        filename = tmp_path / f"{self.test_broken_data_instance.sample_id}_broken_viral_aligner_analysis_fields.json"
        self.test_broken_data_instance.analysis_table.write_analysis_to_json(filename)
        print(f"Saving json to {filename}")
