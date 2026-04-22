"""
utils.py - a collection of functions used by main.py, bacteria.py and virus.py including onyx calls and onyx analysis
table creation function, error handling and io functions.
"""

import json
import logging
import os
from pathlib import Path

import pandas as pd
import yaml
from onyx import OnyxClient, OnyxConfig, OnyxEnv
from onyx_analysis_helper import onyx_analysis_helper_functions as oa
from pandas.core.frame import DataFrame


##############
# Exceptions #
class ClasParError(Exception):
    """Base class for exceptions in this module"""

    pass


class InputError(ClasParError):
    """
    Exception raised for errors in input files.

    :ivar exitcode: int, exitcode
    :ivar error_type: str, error type caught.
    :ivar message: helpful explicit error message.
    """

    def __init__(self, exitcode: int, error_type: str, message: str):
        """
        Custom Error code - supply exitcode (int) and message (str).
        """
        self.exitcode: int = exitcode
        self.message: str = message
        self.error_type: str = error_type

    def __str__(self):
        return repr(self.error_type)


########
# Onyx #

# Set up onyx config
CONFIG = OnyxConfig(
    domain=os.environ[OnyxEnv.DOMAIN],
    token=os.environ[OnyxEnv.TOKEN],
)


@oa.call_to_onyx
def get_input_data(sample_id: str, server: str) -> tuple[list[pd.DataFrame], int]:
    """
    Get the input data from Onyx. Decorated to handle errors suitably.
    :param sample_id: ID of the sample (climb-id).
    :param server: the server to query.
    :return: tuple of dataframes (list) and exitcode (int).
    """
    with OnyxClient(CONFIG) as client:
        record = client.get(
            project=server,
            climb_id=sample_id,
            include=["classifier_calls", "alignment_results", "sylph_results"],
        )

    try:
        alignment_results_df = pd.DataFrame(record["alignment_results"])
        sylph_results_df = pd.DataFrame(record["sylph_results"])
        classifier_calls_df = pd.DataFrame(record["classifier_calls"])
        exitcode = 0
    except KeyError as e:
        logging.error("Could not find key %s in Onyx Record. Exiting cleanly." % (e))  # noqa
        exitcode = 1

    return [alignment_results_df, sylph_results_df, classifier_calls_df], exitcode


###################
# Analysis tables #


def create_analysis_fields(
    *,
    domain: str,
    classifier: str,
    record_id: str,
    thresholds_dict: dict[str, int | str],
    tool_versions: dict,
    headline_result: str,
    results: dict,
    server: str,
) -> tuple[oa.OnyxAnalysis, int]:
    """
    Set up fields dictionary used to populate analysis table containing ClasPar outputs.
    :param domain: str, one of 'bacteria', 'virus', 'fungi' etc
    :param classifier: the type of classifier being reported in the table (kraken or sylph)
    :param record_id: Climb ID for sample
    :param thresholds_dict: Dictionary containing criteria used to filter, which gets added to the
    methods field as 'thresholds': {thresholds_dict}
    :param tool_versions: dict of tools, databases, files etc and their versions.
    :param headline_result: Short description of main result
    :param results: Dictionary containing results
    :param server: Server code is running on, one of "server" or "synthscape"
    :returns: onyx analysis object and exitcode.
    onyx_analysis: Class containing required fields for input to onyx analysis table.
    exitcode: Exit code for checks - will be 0 if all checks passed, 1 if any checks failed
    """
    onyx_analysis = oa.OnyxAnalysis()  # set up class
    # Add analysis details
    onyx_analysis.add_analysis_details(
        analysis_name=f"claspar-{classifier}-{domain}",
        analysis_description=f"This is an analysis to parse and filter the {domain} classifications from {classifier} "
        f"and look up the clinical profiles for classified taxa.",
    )
    # Add metadata about the pipeline/package
    onyx_analysis.add_package_metadata(package_name="claspar")
    # Check that the methods were parsed by the class
    methods_fail = onyx_analysis.add_methods(sample_id=record_id, server_name=server, tool_versions=tool_versions)

    # Reformat the thresholds_dict:
    methods_dict: dict[str, dict[str, int | str]] = {"thresholds": thresholds_dict}
    # Check that additional methods are parsed by the class
    other_methods_fail = onyx_analysis.add_other_methods(methods_dict)

    # Check that the results were parsed by the class
    results_fail = onyx_analysis.add_results(top_result=headline_result, results_dict=results)
    # Add info about sample and server (server/synthscape)
    onyx_analysis.add_server_records(sample_id=record_id, server_name=server)
    # Check the final object using the helper method
    required_field_fail, attribute_fail = onyx_analysis.check_analysis_object(publish_analysis=False)
    # If any fail, raise exit code.
    if any(  # noqa: SIM108
        [methods_fail, other_methods_fail, results_fail, required_field_fail, attribute_fail]
    ):  # noqa SIM108
        exitcode = 1
    else:
        exitcode = 0

    return onyx_analysis, exitcode


#############
# Functions #


def read_samplesheet(path_to_samplesheet: Path | str) -> tuple[list[pd.DataFrame], int]:
    """
    Read in tab seperated samplesheet and return three dataframes. Must be 2x2 dataframe with 'full_Onyx_json' header
    that contains the onyx record as json.

    :param path_to_samplesheet: path to the samplesheet to be read in, must be tab seperated.
    :return: tuple of dataframes (list) and exitcode (int). List of dataframes consists of the alignment results,
    the sylph results and the classifier calls. Note that any of these could be empty dataframes!
    """
    exitcode = 0
    samplesheet_df: DataFrame = pd.read_csv(path_to_samplesheet, sep="\t")
    try:
        json_str = samplesheet_df["full_Onyx_json"].iloc[0]
    except KeyError as k:
        logging.debug("Expected column %s not found in the samplesheet, will try using index..." % (k))  # noqa
        try:
            json_str = samplesheet_df.loc[0, :].values[1]
        except IndexError as i:
            logging.error(
                """
                Cannot parse the json from the sample sheet, dataframe has %s row and %s columns.
                Expected 1 row, 2 columns. Expected columns "climb_id" and "full_Onyx_json". %s
                """  # noqa
                % (samplesheet_df.shape[0], samplesheet_df.shape[1], i)  # noqa
            )

            exitcode = 1
            empty_dfs = [pd.DataFrame() for _ in range(3)]  # Make three empty dataframes to return
            return empty_dfs, exitcode

    # Try to parse the json:
    try:
        record = json.loads(json_str)
    except TypeError as t:
        logging.error("Cannot parse the json from the sample sheet. Expected json, got %s. %s" % (json_str, t))  # noqa
        exitcode = 1
        empty_dfs = [pd.DataFrame() for _ in range(3)]  # Make three empty dataframes to return
        return empty_dfs, exitcode

    alignment_results_df = pd.DataFrame(record["alignment_results"])
    sylph_results_df = pd.DataFrame(record["sylph_results"])
    classifier_calls_df = pd.DataFrame(record["classifier_calls"])

    dfs = [alignment_results_df, sylph_results_df, classifier_calls_df]

    return dfs, exitcode


def check_filters(filters: list[str], threshold_dict: dict) -> int:
    """
    Checks that the filter threshold dictionary contains all the keys as needed. If there are any extra keys, these are
    logged but exitcode is 0. If any keys are missing, these are logged as an error and exitcode 1 is returned.

    :param filters: list of strings: the filters expected in the filter threshold dictionary.
    :param threshold_dict: the filter threshold dictionary.
    :return: exitcode (int)
    """
    exitcode = 0
    # First check for filters in the threshold dict compared to expected (extras)
    for f1 in threshold_dict:
        if f1 not in filters:
            logging.info("Filter '%s' is not set to be used." % (f1))  # noqa

    # Second check for filters not in the threshold dict that are expected (missing)
    for f2 in filters:
        if f2 not in threshold_dict:
            logging.error("Filter %s has not been provided and is required. Exiting." % (f2))  # noqa
            exitcode = 1
    return exitcode


def read_config_file(config_file: str | Path) -> tuple[dict, list]:
    """
    Read config file to get thresholds to filter each of the classifier results. Check that all filters are accounted
    for and log any that are not used or missing.
    :param config_file: path to yaml file containing filter thresholds.
    :returns: dict, nested dictionary of thresholds and list of exit codes (list of ints)
    """
    exit_codes = []

    with Path(config_file).open("r") as file:
        thresholds = yaml.safe_load(file)

    expected_kraken_filters = [
        "READ_THRESHOLD",
        "GENUS_RANK_THRESHOLD",
        "GENUS_READ_PCT_THRESHOLD",
    ]
    exit_codes.append(check_filters(expected_kraken_filters, thresholds["kraken_bacterial_filters"]))

    expected_sylph_filters = [
        "CONTAINMENT_INDEX_THRESHOLD",
        "EFFECTIVE_COVERAGE_THRESHOLD",
    ]
    exit_codes.append(check_filters(expected_sylph_filters, thresholds["sylph_filters"]))

    expected_viral_aligner_filters = [
        "EVENNESS_VALUE",
        "COVERAGE_1X",
        "UNIQUELY_MAPPED_READS",
        "MEAN_READ_IDENTITY",
        "MEAN_ALIGNMENT_LENGTH",
    ]
    exit_codes.append(check_filters(expected_viral_aligner_filters, thresholds["viral_aligner_filters"]))
    return thresholds, exit_codes


def setup_outdir(outdir: str | Path) -> None:
    """
    Handle the output directory, using either the default or commandline arg. Create dir
    if needed.
    :param outdir: the outdir argument from commandline.
    :return: None.
    """
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)

    return None


def write_df_to_csv(*, df: pd.DataFrame, filename: str, results_dir: str | Path) -> Path:
    """
    Write results dataframe to csv.
    :param df: dataframe of results to save.
    :param filename: str, unique name of file WITHOUT extension. (csv gets added).
    :param results_dir: Directory to save results to.
    :returns: os.path of saved csv file.
    """

    result_file_path = Path(results_dir) / f"{filename}.csv"

    df.to_csv(result_file_path, index=False)

    return result_file_path
