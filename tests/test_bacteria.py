"""
Create unit tests for modules in the tests/ folder. All functions in a repo should be unit tested and
tests should be run before and after any changes are made.
"""

import json

import pandas as pd
import pytest  # noqa: F401
from taxaplease import TaxaPlease

from claspar import bacteria

pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)


tp = TaxaPlease()


@pytest.fixture(scope="module")
def metadata_json():
    with open("tests/test_unit_test_metadata.json") as json_file:  # noqa: PTH123
        metadata = json.loads(json_file.read())
    assert metadata["sylph_results"] and metadata["classifier_calls"]
    return metadata


def test_input(metadata_json):
    print(metadata_json["classifier_calls"])


################
# Kraken Tests #
################


class TestKrakenBacteria:
    @pytest.fixture(autouse=True)
    def kraken_thresholds_dict(self):
        self.thresholds = {
            "READ_THRESHOLD": 10,
            "GENUS_RANK_THRESHOLD": 3,
            "GENUS_READ_PCT_THRESHOLD": 20,
        }

    @pytest.fixture(autouse=True)
    def test_kraken_data(self, metadata_json):
        self.test_input_data = pd.DataFrame(metadata_json["classifier_calls"])

    @pytest.fixture(autouse=True)
    def test_instance_1(self, kraken_thresholds_dict, test_kraken_data):
        kraken_assignments = bacteria.KrakenBacteria(
            sample_id="ID-12345678",
            original_classifier_df=self.test_input_data,
            kraken_bacteria_thresholds_dict=self.thresholds,  # ty:ignore[invalid-argument-type]
            profiles_dict={
                "ProfileX": {
                    662: {"taxon": "Vibrio", "rank": "genus"},
                    547: {"taxon": "Enterobacteriaceae", "rank": "family"},
                },
                "ProfileY": {
                    22: {"taxon": "Shewanella", "rank": "genus"},
                    520: {"taxon": "Bordetella pertussis", "rank": "species"},
                },
                "ProfileZ": {
                    517: {"taxon": "Bordetella", "rank": "genus"},
                },
            },
            taxaplease_instance=tp,
            server="server",
        )
        self.kraken_class_instance = kraken_assignments

    def test_test_class(self):
        assert (s := self.kraken_class_instance.sample_id) == "ID-12345678", (
            f"Sanity checking the test class setup failed - expected ID-12345678, got {s}"
        )

    @pytest.mark.parametrize(
        "taxonid,isbacteria,parent_id,parent_name",
        [
            (139, True, 64895, "Borreliella"),
            (1410656, True, 859, "Fusobacterium necrophorum"),
            (2696357, False, 2788787, "unclassified Caudoviricetes"),
            (3052230, False, 11102, "Hepacivirus"),
        ],
    )
    def test__get_parent_taxonomy(self, taxonid, isbacteria, parent_id, parent_name):
        actual = self.kraken_class_instance._get_parent_taxonomy(taxonid)
        assert (a := actual[0]) == isbacteria, f"Expected {isbacteria} but got {a}"
        assert (b := actual[1]["taxid"]) == parent_id, f"Expected {parent_id} but got {b}"  # type: ignore
        assert (c := actual[2]) == parent_name, f"Expected {parent_name}, but got {c}"
        print(f"\nParent taxonomy dict for {taxonid} is {actual[1]}, all is as expected.")

    @pytest.mark.parametrize(
        "test,count_descendants,order_in_genus,pct_genus_reads,outcome",
        [
            ("Dominant taxa, high count", 100, 1, 80, "high"),
            ("thresholds", 10, 3, 20, "high"),
            ("low count, many species in genus", 5, 5, 5, "low"),
            ("high count, many species in genus", 1000, 1, 10, "low"),
        ],
    )
    def test__get_kraken_confidence_rating(self, test, count_descendants, order_in_genus, pct_genus_reads, outcome):
        actual = self.kraken_class_instance._get_kraken_confidence_rating(
            count_descendants=count_descendants, order_in_genus=order_in_genus, pct_genus_reads=pct_genus_reads
        )
        assert actual == outcome, f"Expected {outcome}, got {actual}"
        print(
            f"\nTesting '{test}' with {count_descendants} reads, {order_in_genus} rank in the genus, "
            f"{pct_genus_reads} percentage genus reads was given {outcome} confidence rating. All is as expected."
        )

    def test__add_pct_and_rank(self):
        species_df = pd.DataFrame(
            [
                [294, "Pseudomonas fluorescens", 0, 15, 15, 286, 100, 15, 2],
                [317, "Pseudomonas syringae", 0, 15, 15, 286, 100, 15, 2],
                [28450, "Burkholderia pseudomallei", 0, 15, 14, 32008, 50, 30, 2],
                [33069, "Pseudomonas viridiflava", 0, 15, 15, 286, 100, 15, 2],
                [60550, "Burkholderia pyrrocinia", 0, 20, 20, 32008, 50, 40, 1],
                [95485, "Burkholderia stabilis", 0, 15, 14, 32008, 50, 30, 2],
                [190893, "Vibrio coralliilyticus", 0.87, 50, 50, 662, 200, 25, 2],
                [190897, "Vibrio gallicus", 1.05, 60, 60, 662, 200, 30, 1],
                [553611, "Photobacterium leiognathi", 1.01, 500, 0, 657, 500, 100, 1],
                [1637837, "Burkholderia savannae", 0, 2, 2, 32008, 50, 4, 3],
                [2163016, "Vibrio sp. Dhg", 0.45, 40, 40, 662, 200, 20, 3],
                [2200953, "Vibrio albus", 0.22, 20, 2, 662, 200, 10, 4],
                [2496559, "Vibrio aquaticus", 0.12, 10, 5, 662, 200, 5, 5],
                [2961893, "Pseudomonas sp. ZM23", 0, 15, 15, 286, 100, 15, 2],
                [3042029, "Pseudomonas sp. PMCC200367", 0, 40, 40, 286, 100, 40, 1],
            ],
            columns=[
                "taxon_id",
                "human_readable",
                "percentage",
                "count_descendants",
                "count_direct",
                "genus_id",
                "genus_level_reads",
                "expected_pct_genus",
                "expected_rank",
            ],
        )
        species_df = self.kraken_class_instance._add_pct_and_rank(species_df)
        assert all(species_df["pct_genus_reads"] == species_df["expected_pct_genus"])
        assert all(species_df["order_in_genus"] == species_df["expected_rank"])

    def test__process_kraken(self):
        actual_species_df, actual_genus_df = self.kraken_class_instance._process_kraken()
        total_species = len(actual_species_df)

        assert total_species == 12, f"Expected 12, got {total_species}"

        high_confidence_species = actual_species_df.loc[actual_species_df["kraken_confidence"] == "high"].shape[0]
        assert high_confidence_species == 8, f"Expected 8, got {high_confidence_species}"

        total_genera = len(actual_genus_df)
        assert total_genera == 7, f"Expected 7, got {total_genera}"

        genera_with_more_than_one_species = actual_genus_df.loc[actual_genus_df["total_species_identified"] > 1].shape[
            0
        ]
        assert genera_with_more_than_one_species == 3, f"Expected 3, got {genera_with_more_than_one_species}"

        print(
            f"\nProcessing kraken results for the test data revealed {total_species} total species, of which "
            f"{high_confidence_species} were high confidence. There were {total_genera} total genera, of which "
            f"{genera_with_more_than_one_species} had more than one species in it. (All as expected)"
        )

    def test_get_kraken_results(self):
        headline_result, main_result, species_df, genus_df = self.kraken_class_instance._get_kraken_results()
        expected_headline = (
            "Kraken classified 12 bacterial species and 2 genera; ProfileX-high, ProfileY-high, ProfileY-low, "
            "Unknown-high, Unknown-low"
        )

        assert expected_headline in headline_result
        assert all(len(inner.keys()) == 4 for inner in main_result.values())

        assert len(main_result) == 4, f"Expected 4 taxa in the final result, got {len(main_result)}"

        print(headline_result)
        print(main_result)

    def test_genera_handled_correctly(self):
        headline_result, main_result, species_df, genus_df = self.kraken_class_instance._get_kraken_results()
        # these genera have no species classified
        expected_genera_in_results = ["Nostoc", "Calothrix"]

        assert all(g for g in expected_genera_in_results if g in genus_df["human_readable"])

    def test_save_outputs_to_csv(self, tmp_path):
        self.kraken_class_instance.save_outputs_to_csv(tmp_path)
        print(f"Saving to {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.kraken_class_instance.get_kraken_bacteria_analysis_table("test_profile_table.xlsx")

        filename = tmp_path / f"{self.kraken_class_instance.sample_id}_kraken_bacteria_analysis_fields.json"
        self.kraken_class_instance.analysis_table.write_analysis_to_json(filename)
        print(f"Saving json to {filename}")


class TestNoKrakenBacteria:
    @pytest.fixture(autouse=True)
    def kraken_thresholds_dict(self):
        self.thresholds: dict[str, int] = {
            "READ_THRESHOLD": 10,
            "GENUS_RANK_THRESHOLD": 3,
            "GENUS_READ_PCT_THRESHOLD": 20,
        }

    @pytest.fixture(autouse=True)
    def test_kraken_data(self):
        self.test_input_data = pd.DataFrame()

    @pytest.fixture(autouse=True)
    def test_instance_no_data(self, kraken_thresholds_dict, test_kraken_data):
        kraken_assignments = bacteria.KrakenBacteria(
            sample_id="ID-00000000",
            original_classifier_df=self.test_input_data,
            kraken_bacteria_thresholds_dict=self.thresholds,  # ty:ignore[invalid-argument-type]
            profiles_dict={},
            taxaplease_instance=tp,
            server="server",
        )
        self.kraken_class_instance = kraken_assignments

    def test_test_class(self):
        assert (s := self.kraken_class_instance.sample_id) == "ID-00000000", (
            f"Sanity checking the test class setup failed - expected ID-0000000, got {s}"
        )

    def test__process_kraken(self):
        species_df, genus_df = self.kraken_class_instance._process_kraken()
        assert species_df.empty, f"Expected empty species df, got {species_df}"
        assert genus_df.empty, f"Expected empty genus df, got {genus_df}"

    def test__get_kraken_results(self):
        headline, result, species_df, genus_df = self.kraken_class_instance._get_kraken_results()
        assert "Kraken classified 0 bacterial species." in headline
        assert result == {}, f"Expected empty dict, got {result}"
        assert species_df.empty, f"Expected empty species df, got {species_df}"
        assert genus_df.empty, f"Expected empty genus df, got {genus_df}"

    def test_save_outputs_to_csv(self, tmp_path):
        self.kraken_class_instance.save_outputs_to_csv(tmp_path)
        print(f"Saving to {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.kraken_class_instance.get_kraken_bacteria_analysis_table("test_profile_table.xlsx")

        filename = tmp_path / f"{self.kraken_class_instance.sample_id}_kraken_bacteria_analysis_fields.json"
        self.kraken_class_instance.analysis_table.write_analysis_to_json(filename)
        print(f"Saving json to {filename}")


###############
# Sylph Tests #
###############


class TestSylphBacteria:
    @pytest.fixture(autouse=True)
    def sylph_test_data(self, metadata_json):
        self.sylph_test_df = pd.DataFrame(metadata_json["sylph_results"])

    @pytest.fixture(autouse=True)
    def sylph_thresholds_dict(self):
        self.thresholds = {"CONTAINMENT_INDEX_THRESHOLD": 0.2, "EFFECTIVE_COVERAGE_THRESHOLD": 1.0}

    @pytest.fixture(autouse=True)
    def sylph_test_data_with_rank(self, sylph_test_data):
        sylph_test_data = self.sylph_test_df.copy()
        sylph_test_data["taxon_rank"] = sylph_test_data["taxon_id"].apply(lambda x: tp.get_record(x)["rank"])  # type: ignore
        self.sylph_df_with_rank = sylph_test_data

    @pytest.fixture(autouse=True)
    def test_instance_1(self, sylph_test_data, sylph_thresholds_dict):
        sylph_assignments_1 = bacteria.SylphBacteria(
            sample_id="ID-12345678",
            original_sylph_df=self.sylph_test_df,
            sylph_bacteria_thresholds_dict=self.thresholds,
            profiles_dict={
                "ProfileX": {
                    520: {"taxon": "Bordetella pertussis", "rank": "species"},
                },
                "ProfileY": {
                    2104: {"taxon": "Mycoplasmoides pneumoniae", "rank": "species"},
                    520: {"taxon": "Bordetella pertussis", "rank": "species"},
                },
                "ProfileZ": {
                    517: {"taxon": "Bordetella", "rank": "genus"},
                },
            },
            taxaplease_instance=tp,
            server="server",
        )
        self.sylph_class_instance_1 = sylph_assignments_1

    def test_test_class(self):
        assert (s := self.sylph_class_instance_1.sample_id) == "ID-12345678", (
            f"Sanity checking the test class setup failed - expected ID-12345678, got {s}"
        )

    def test__process_sylph_rank(self):
        expected_new_columns = pd.DataFrame(
            [[520, "Bordetella pertussis"], [2104, "Mycoplasmoides pneumoniae"], [519, "Bordetella parapertussis"]],
        )

        actual_new_columns = self.sylph_df_with_rank.apply(
            lambda x: self.sylph_class_instance_1._process_sylph_rank(x), axis=1, result_type="expand"
        )
        assert actual_new_columns.equals(expected_new_columns)
        print(f"\nChecking the sylph rank to get ID and species name - get {actual_new_columns}")

    def test__process_sylph_rank_strain(self):
        test_df = pd.DataFrame(
            {"taxon_id": [1121296], "human_readable": ["[Clostridium] aminophilum DSM 10710"], "taxon_rank": ["strain"]}
        )
        expected_new_columns = pd.DataFrame([[1526, "[Clostridium] aminophilum"]])

        actual_new_columns = test_df.apply(
            lambda x: self.sylph_class_instance_1._process_sylph_rank(x), axis=1, result_type="expand"
        )
        assert actual_new_columns.equals(expected_new_columns)
        print(
            f"\nTaxa with strain rank will be given the species level ID and human readable:"
            f"{test_df} would give {expected_new_columns}."
        )

    def test__process_sylph_rank_no_id(self):
        test_df = pd.DataFrame({"taxon_id": [""], "human_readable": [""], "taxon_rank": [""]})
        expected_new_columns = pd.DataFrame([[None, None]])

        actual_new_columns = test_df.apply(
            lambda x: self.sylph_class_instance_1._process_sylph_rank(x), axis=1, result_type="expand"
        )
        assert actual_new_columns.equals(expected_new_columns)
        print(
            f"\nTaxa without a rank will be given the species level ID and human readable:"
            f"{test_df} would give {expected_new_columns}."
        )

    @pytest.mark.parametrize(
        "test,containment_index,effective_coverage,outcome",
        [
            ("high containment, high coverage", 1.0, 10.0, "high"),
            ("thresholds", 0.2, 1.0, "high"),
            ("low containment, high coverage", 0.05, 0.5, "low"),
            ("low containment, low coverage", 0.1, 0.04, "low"),
        ],
    )
    def test__get_sylph_confidence_rating(self, test, containment_index, effective_coverage, outcome):
        actual_outcome = self.sylph_class_instance_1._get_sylph_confidence_rating(
            containment_index=containment_index, effective_coverage=effective_coverage
        )

        assert actual_outcome == outcome
        print(
            f"\nTest '{test}' with containment index {containment_index} and effective coverage {effective_coverage} "
            f"has confidence rating {outcome} (as expected)."
        )

    def test__process_sylph(self):
        actual_sylph_result = self.sylph_class_instance_1._process_sylph()
        print(actual_sylph_result)
        species_count = actual_sylph_result.loc[actual_sylph_result["taxon_rank"] == "species"].shape[0]
        assert species_count == 2  # 2 species level assignments, 1 strain
        high_confidence_sylph = actual_sylph_result.loc[actual_sylph_result["sylph_confidence"] == "high"].shape[0]
        assert high_confidence_sylph == 2  # 2 high confidence level taxa.
        count_of_rows_with_genus_ids = actual_sylph_result.loc[actual_sylph_result["genus_id"].notnull()].shape[0]
        # Genus ID is null if there is no genus for the species, for example if it's unclassified.
        assert count_of_rows_with_genus_ids == 3
        print(
            f"\nProcessing sylph data - there were {len(actual_sylph_result)} sylph results, of which {species_count}"
            f" were species, and {high_confidence_sylph} had a high confidence. There were "
            f"{count_of_rows_with_genus_ids} classifications with a genus ID (all as expected)."
        )

    def test__get_sylph_results_instance1(self):
        headline_result, results, sylph_processed_df = self.sylph_class_instance_1._get_sylph_results()
        print(headline_result)
        expected_headline = "Sylph classified 3 bacterial (or archael) taxa; ProfileX-high, ProfileY-high, ProfileZ-low"
        assert expected_headline in headline_result

        assert len(results.keys()) == 2  # 2 taxa with high confidence.
        assert all(len(inner.keys()) == 4 for inner in results.values())
        assert sylph_processed_df.shape == (3, 26)  # 3 processed taxa

    def test_get_sylph_analysis_table(self):
        analysis_table = self.sylph_class_instance_1.get_sylph_analysis_table("test_profile_table.xlsx")
        assert (p := analysis_table.pipeline_name) == "ClasPar", f'Expected pipeline name "ClasPar", got "{p}"'
        assert (n := analysis_table.name) == "claspar-sylph-bacteria", (
            f'Expected name "claspar-sylph-bacteria", got {n}'
        )
        assert "sylph" in (d := analysis_table.description), f'Expected "sylph" to be in the description, got {d}'

    def test_save_outputs_to_csv(self, tmp_path):
        self.sylph_class_instance_1.save_outputs_to_csv(tmp_path)
        print(f"Saving to {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.sylph_class_instance_1.get_sylph_analysis_table("test_profile_table.xlsx")

        filename = tmp_path / f"{self.sylph_class_instance_1.sample_id}_sylph_analysis_fields.json"
        self.sylph_class_instance_1.analysis_table.write_analysis_to_json(filename)
        print(f"Saving json to {filename}")

    ###### Make another instance of the class with some different data:
    @pytest.fixture(autouse=True)
    def test_instance_2(self, sylph_test_data, sylph_thresholds_dict):
        sylph_test_df_edited = self.sylph_test_df.copy()
        # add a match at the genus level
        sylph_test_df_edited.loc[3] = [
            "838",
            "Prevotella",
            "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella ruminicola_D",
            "GCA_946637185.1",
            "SRR873610_bin.90_metaWRAP_v1.3_MAG",
            71.22,
            61.552,
            98.14,
            "NA-NA",
            17.451,
            "NA-NA",
            15,
            8.991,
            "2111/2641",
            98.14,
            0,
        ]
        # add a match at the family level:
        sylph_test_df_edited.loc[4] = [
            "186803",
            "Lachnospiraceae",
            "d__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Butyrivibrio;s__Butyrivibrio fibrisolvens",
            "GCF_037113525.1",
            "CP146963.1 Butyrivibrio fibrisolvens strain D1 chromosome, complete genome",
            48.65,
            47.125,
            96.66,
            "NA-NA",
            14.879,
            "NA-NA",
            14,
            7.619,
            "87/1992",
            96.79,
            "0",
        ]
        self.sylph_test_df = sylph_test_df_edited

        sylph_assignments_2 = bacteria.SylphBacteria(
            sample_id="ID-87654321",
            original_sylph_df=self.sylph_test_df,
            sylph_bacteria_thresholds_dict=self.thresholds,
            profiles_dict={
                "ProfileX": {
                    520: {"taxon": "Bordetella pertussis", "rank": "species"},
                    838: {"taxon": "Prevotella", "rank": "genus"},
                },
                "ProfileY": {
                    2104: {"taxon": "Mycoplasmoides pneumoniae", "rank": "species"},
                    520: {"taxon": "Bordetella pertussis", "rank": "species"},
                    186803: {"taxon": "Lachnospiraceae", "rank": "family"},
                },
                "ProfileZ": {
                    517: {"taxon": "Bordetella", "rank": "genus"},
                },
            },
            taxaplease_instance=tp,
            server="server",
        )
        self.sylph_class_instance_2 = sylph_assignments_2

    def test_test_class_again(self):
        assert (s := self.sylph_class_instance_2.sample_id) == "ID-87654321", (
            f"Sanity checking the test class setup failed - expected ID-87654321, got {s}"
        )

    def test__process_sylph_with_genus(self):
        actual_sylph_result = self.sylph_class_instance_2._process_sylph()
        # The penultimate row (index = 3) is Prevotella at genus level.
        assert pd.isna(actual_sylph_result.at[3, "species_id"])
        assert actual_sylph_result.at[3, "taxon_rank"] == "genus"  # check the 4th row in the outputs.
        assert actual_sylph_result.at[3, "genus_id"] == "838"
        # The last row (index = 4) is Butyrivibrio at the family (Lachnospiraceae) level.
        assert pd.isna(actual_sylph_result.at[4, "species_id"])
        assert actual_sylph_result.at[4, "taxon_rank"] == "family"
        assert pd.isna(actual_sylph_result.at[4, "genus_id"])

        print(actual_sylph_result)
        assert not pd.isna(actual_sylph_result.at[3, "genus_id"]), f"{actual_sylph_result.at[3, 'genus_id']}"

    def test__get_sylph_results_instance2(self):
        headline_result, results, sylph_processed_df = self.sylph_class_instance_2._get_sylph_results()
        # print(sylph_processed_df)
        # print(headline_result)
        expected_headline = "Sylph classified 5 bacterial (or archael) taxa; ProfileX-high, ProfileY-high, ProfileY-low"
        assert expected_headline in headline_result
        # print(results)
        assert len(results.keys()) == 3
        assert all(len(inner.keys()) == 4 for inner in results.values())

        assert sylph_processed_df.shape[0] == 5


class TestNoSylphBacteria:
    @pytest.fixture(autouse=True)
    def sylph_test_data(self):
        self.sylph_test_df = pd.DataFrame()

    @pytest.fixture(autouse=True)
    def sylph_thresholds_dict(self):
        self.thresholds = {"CONTAINMENT_INDEX_THRESHOLD": 0.2, "EFFECTIVE_COVERAGE_THRESHOLD": 1.0}

    @pytest.fixture(autouse=True)
    def test_instance_1(self, sylph_test_data, sylph_thresholds_dict):
        sylph_assignments = bacteria.SylphBacteria(
            sample_id="ID-00000000",
            original_sylph_df=self.sylph_test_df,
            sylph_bacteria_thresholds_dict=self.thresholds,
            profiles_dict={},
            taxaplease_instance=tp,
            server="server",
        )
        self.sylph_class_instance = sylph_assignments

    def test_test_class(self):
        assert (s := self.sylph_class_instance.sample_id) == "ID-00000000", (
            f"Sanity checking the test class setup failed - expected ID-00000000, got {s}"
        )

    def test__get_sylph_results(self):
        headline, result, processed_sylph_df = self.sylph_class_instance._get_sylph_results()
        assert "Sylph classified 0 bacterial (or archaeal) taxa." in headline
        assert result == {}, f"Expected empty dict, got {result}"
        assert processed_sylph_df.empty, f"Expected empty genus df, got {processed_sylph_df}"

    def test_process_sylph(self):
        actual_sylph_result = self.sylph_class_instance._process_sylph()
        assert actual_sylph_result.empty, f"Expected empty dataframe, got {actual_sylph_result}"

    def test_save_outputs_to_csv(self, tmp_path):
        self.sylph_class_instance.save_outputs_to_csv(tmp_path)
        print(f"Saving to {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.sylph_class_instance.get_sylph_analysis_table("test_profile_table.xlsx")

        filename = tmp_path / f"{self.sylph_class_instance.sample_id}_sylph_analysis_fields.json"
        self.sylph_class_instance.analysis_table.write_analysis_to_json(filename)
        print(f"Saving json to {filename}")
