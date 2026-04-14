"""
Module for handling the viral aligner outputs. Filters are applied to return the taxa of interest.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from onyx_analysis_helper import onyx_analysis_helper_functions as oa
from profiler import get_profiles as profiler
from taxaplease import TaxaPlease

from claspar.utils import create_analysis_fields, write_df_to_csv


class VirusClasPar:
    """
    Class for parsing the viral aligner classifications.

    The main methods are:
    get_virus_analysis_table - get the analysis table (returns instance of the onyx analysis helper class).
    save_outputs_to_csv - method to save the outputs (saved in instance attributes) to file (returns None)

    Attributes:
    data_input - pd.DataFrame; original viral aligner results from scylla.
    thresholds - dict; the thresholds for filtering
    sample_id - str; the climb-id
    server - str;  serveror synthscape
    profiles_dict - dict; lookup for profiles and their taxa.
    taxaplease - instance of taxaplease. Will instantiate a new taxaplease instance one if not given.
    filtered_data - pd.DataFrame; the viral aligner results after filtering
    headline_results - str; the main result, automatically generated to include the final number of taxa that remained
    after filtering
    results - dict; the filtered dataframe as a dict.
    analysis_table - oa.OnyxAnalysis; instance of the analysis table from the helper, containing all the relevant info.

    :param sample_id: str, climb-id
    :param original_viral_aligner_df: pandas dataframe, the original results from scylla.
    :param virus_thresholds_dict: dict, containing the thresholds to filter.
    :param profiles_dict: dict containing profiles and their taxa.
    :param taxaplease_instance: instance of TaxaPlease class (optional).
    :param server: str,  serveror synthscape (database server for Onyx to connect to - will be validated.)
    """

    def __init__(
        self,
        sample_id: str,
        original_viral_aligner_df: pd.DataFrame,
        virus_thresholds_dict: dict[str, int | float],
        profiles_dict: dict,
        taxaplease_instance: TaxaPlease | None,
        server: str = "server",
    ):
        """
        Create instance of VirusClasPar class, where arguments are attributes and instance methods populate
        headline_result, results and sylph_processed_df.

        :param sample_id: str, climb-id
        :param original_classifier_df: pandas dataframe, the original results from scylla.
        :param kraken_bacteria_thresholds_dict: dict, containing the thresholds to filter.
        :param profiles_dict: dict containing profiles and their taxa.
        :param taxaplease_instance: instance of TaxaPlease class (optional).
        :param server: str, database server for Onyx to connect to - will be validated.
        """
        self.viral_aligner_df: pd.DataFrame = original_viral_aligner_df
        self.thresholds: dict = virus_thresholds_dict
        self.sample_id: str = sample_id
        self.profiles_dict: dict[str, dict[int, dict[str, str]]] = profiles_dict
        self.taxaplease: TaxaPlease = taxaplease_instance if taxaplease_instance else TaxaPlease()
        self.server: str = server

        self.headline_results: str
        self.results: dict
        self.processed_viral_aligner_df: pd.DataFrame
        self.headline_results, self.results, self.processed_viral_aligner_df = self._get_viral_aligner_results()

    def _process_viral_aligner(self) -> pd.DataFrame:
        """
        Process the viral aligner data by applying the filters to assign low or high confidence.

        :param viral_aligner_df: dataframe of the viral aligner outputs straight from scylla.
        :param thresholds_dict: dictionary of the filters to apply to the dataframe.
        :return: filtered df. Filtered df is pandas dataframe with the high confidence taxa.
        (Empty df returned if input is also empty.)
        """
        if self.viral_aligner_df.empty:
            return pd.DataFrame()

        processed_viral_aligner_df = self.viral_aligner_df.copy()  # Don't edit the original dataframe

        # Add in the taxon_rank:
        processed_viral_aligner_df["taxon_rank"] = processed_viral_aligner_df["taxon_id"].apply(
            lambda x: rec["rank"] if (rec := self.taxaplease.get_record(x)) is not None else "Unknown"
        )

        try:
            processed_viral_aligner_df["confidence"] = np.where(
                (processed_viral_aligner_df["evenness_value"] >= self.thresholds["EVENNESS_VALUE"])
                & (processed_viral_aligner_df["coverage_1x"] >= self.thresholds["COVERAGE_1X"])
                & (processed_viral_aligner_df["uniquely_mapped_reads"] >= self.thresholds["UNIQUELY_MAPPED_READS"])
                & (processed_viral_aligner_df["mean_read_identity"] >= self.thresholds["MEAN_READ_IDENTITY"])
                & (processed_viral_aligner_df["mean_alignment_length"] >= self.thresholds["MEAN_ALIGNMENT_LENGTH"]),
                "high",
                "low",
            )
            return processed_viral_aligner_df

        except KeyError as k:
            msg = (
                "WARNING: The viral alignment data from Scylla is missing expected column:%s."
                "Skipping Viral Aligner." % (k)
            )
            logging.error(msg)
            return pd.DataFrame()

    def _get_viral_aligner_results(self) -> tuple[str, dict, pd.DataFrame]:
        """
        Get the main viral aligner results, add profiles, get the headline result and the results-metrics from viral
        aligner outputs.

        :return: tuple; headline_result (str), result (dict), processed_df (dataframe).
        If no alignment data, headline_result is string and results is empty dict.
        """
        processed_va_df = self._process_viral_aligner()

        if processed_va_df.empty:
            headline_result = "Viral Aligner classified 0 viral taxa."
            result = {}
            return headline_result, result, processed_va_df

        # Add profiles
        processed_va_df = profiler.add_profile_to_results(processed_va_df, self.profiles_dict, self.taxaplease)

        # Make a table of the unique profiles with high confidence for the results:
        unique_profile_taxa = processed_va_df.sort_values(["profile", "confidence"])
        unique_profile_taxa = unique_profile_taxa[
            ["profile_taxon_match", "profile_taxon_id", "profile", "confidence"]
        ].drop_duplicates()

        # Results is a dict of the profile taxa that were matched and only high confidence:
        results = (
            unique_profile_taxa.loc[unique_profile_taxa["confidence"] == "high"]
            .dropna()
            .reset_index(drop=True)
            .to_dict(orient="index")
        )
        #### Make headline result
        # Make headline result with format "ProfileA-high, ProfileA-low".

        profile_confidence = unique_profile_taxa.apply(lambda row: f"{row['profile']}-{row['confidence']}", axis=1)
        profile_confidence = sorted(set(profile_confidence.to_list()))
        total_va_taxa = processed_va_df.shape[0]

        headline_result = f"Viral Aligner classified {total_va_taxa} viral taxa; {', '.join(profile_confidence)}"

        return headline_result, results, processed_va_df

    def get_virus_analysis_table(self, profile_table_name: str) -> oa.OnyxAnalysis:
        """
        Pull together all the class attributes into the analysis table.

        :param profile_table_name: Name of the profile table used; this acts as a version.
        """

        analysis_table, _ = create_analysis_fields(
            domain="virus",
            classifier="viralaligner",
            record_id=self.sample_id,
            thresholds=self.thresholds,
            profile_table_name=profile_table_name,
            headline_result=self.headline_results,
            results=self.results,
            server=self.server,
        )

        # Populate the instance attribute:
        self.analysis_table: oa.OnyxAnalysis = analysis_table

        return analysis_table

    def save_outputs_to_csv(self, results_dir: str | Path) -> None:
        """
        Save the final results to csv.
        :param filename: str, name of file to save to.
        :param results_dir: str or path to directory to save to.
        """
        unique_filename = f"{self.sample_id}_viral_aligner_processed"
        write_df_to_csv(df=self.processed_viral_aligner_df, filename=unique_filename, results_dir=results_dir)
