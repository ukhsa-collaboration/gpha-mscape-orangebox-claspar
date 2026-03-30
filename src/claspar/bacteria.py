"""
Module containing functions needed to parse, filter and create bacteria classifications analysis table.
"""

###########
# Imports #
import logging
from pathlib import Path

import pandas as pd
from profiler import get_profiles as profiler  # ty:ignore[unresolved-import]
from taxaplease import TaxaPlease

from claspar.utils import create_analysis_fields, write_df_to_csv


##########
# Kraken #
##########
class KrakenBacteria:
    """
    Class for parsing the bacterial classifications from Kraken.

    The main methods are:
    get_kraken_bacteria_analysis_table - get the analysis table (returns instance of the onyx analysis helper class).
    save_outputs_to_csv - method to save the outputs (saved in instance attributes) to file (returns None)

    Atrributes:
    classifier_results - the original classifier outputs from Scylla.
    thresholds - dict; the thresholds to filter on.
    sample_id - str; climb-id
    profiles_dict - dict; lookup for profiles and their taxa.
    taxaplease - instance of taxaplease. Will instantiate a new taxaplease instance one if not given.
    server - str; databse server.
    kraken_species_results - pd.DataFrame: All the species kraken identified for the sample, plus the genus id and reads
     at genus level, total species in genus identified (and species that pass the filters), the proportion of total
     genus reads, the rank of that in its genus and the kraken confidence (high or low).
    kraken_genus_results - pd.DataFrame: All the genera kraken identified for the sample, plus some info about the
     species within the genus.
    headline_results - str; the main result, automatically generated to include the final number of taxa that were
     assigned high confidence.
    results - dict; the kraken_species_results dataframe filtered to high confidence species as a dict.
    analysis_table - oa.OnyxAnalysis; instance of the analysis table from the helper, containing all the relevant info.


    :param sample_id: str, climb-id
    :param original_classifier_df: pandas dataframe, the original results from scylla.
    :param kraken_bacteria_thresholds_dict: dict, containing the thresholds to filter.
    :param profiles_dict: dict containing profiles and their taxa.
    :param taxaplease_instance: instance of TaxaPlease class (optional).
    :param server: str, database server for Onyx to connect to - will be validated.
    """

    def __init__(
        self,
        sample_id: str,
        original_classifier_df: pd.DataFrame,
        kraken_bacteria_thresholds_dict: dict[str, int | float],
        profiles_dict: dict,
        taxaplease_instance: TaxaPlease | None,
        server: str = "server",
    ):
        """
        Create instance of KrakenBacteria class, where arguments are attributes and instance methods populate
        headline_result, results and kraken_species_results and kraken_genus_results.

        :param sample_id: str, climb-id
        :param original_classifier_df: pandas dataframe, the original results from scylla.
        :param kraken_bacteria_thresholds_dict: dict, containing the thresholds to filter.
        :param profiles_dict: dict containing profiles and their taxa.
        :param taxaplease_instance: instance of TaxaPlease class (optional).
        :param server: str, database server for Onyx to connect to - will be validated.
        """

        self.classifier_results: pd.DataFrame = original_classifier_df
        self.thresholds: dict = kraken_bacteria_thresholds_dict
        self.sample_id: str = sample_id
        self.profiles_dict: dict[str, dict[int, dict[str, str]]] = profiles_dict
        self.taxaplease: TaxaPlease = taxaplease_instance if taxaplease_instance else TaxaPlease()
        self.server: str = server

        self.headline_result: str
        self.result: dict

        self.kraken_species_results: pd.DataFrame
        self.kraken_genus_results: pd.DataFrame

        self.headline_result, self.results, self.kraken_species_results, self.kraken_genus_results = (
            self._get_kraken_results()
        )

    def _get_parent_taxonomy(self, taxon_id: int) -> tuple[bool, dict | None, str | None]:
        """
        Use taxaplease to get 'is bacteria', 'parent record' and 'parent human-readable'.

        :param taxon_id: int, taxon ID
        :return: boolean ('is bacteria'), dict ('parent record') and str ('parent human-readable').
        """
        tp = self.taxaplease

        if taxon_id == 0:
            is_bacteria = False
            parent_record = {}
            parent_human_readable = ""
            return is_bacteria, parent_record, parent_human_readable

        else:
            is_bacteria = tp.isBacteria(taxon_id)
            parent_record = tp.get_parent_record(taxon_id)
            parent_human_readable = parent_record.get("name", None) if parent_record else ""
            return is_bacteria, parent_record, parent_human_readable

    def _get_kraken_confidence_rating(self, *, count_descendants, order_in_genus, pct_genus_reads):
        """
        Apply the thresholds to determine the kraken confidence rating.
        """
        if (
            count_descendants >= self.thresholds["READ_THRESHOLD"]
            and order_in_genus <= self.thresholds["GENUS_RANK_THRESHOLD"]
            and pct_genus_reads >= self.thresholds["GENUS_READ_PCT_THRESHOLD"]
        ):
            return "high"
        else:
            return "low"

    @staticmethod
    def _add_pct_and_rank(species_df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate the proportion of the species reads from the total of genus reads (as percentage).
        Then order all species in a genus but the proportion and assign the rank. If multiple species have the same
        percentage, they get the same rank (the higher number).

        Returns modified dataframe with added columns "pct_genus_reads" and "order_in_genus".

        :param df: dataframe - species dataframe, must have columns "count_descendants" and "genus_level_reads".
        :type df: pd.DataFrame

        """
        df = species_df.copy()
        # calculate percentage of genus level reads
        df["pct_genus_reads"] = (
            df["count_descendants"] / df["genus_level_reads"] * 100
        )  # Proportion of the species reads from total genus reads
        df["order_in_genus"] = (
            df.groupby("genus_id")["pct_genus_reads"].rank(method="dense", ascending=False).astype(int)
        )
        # order by the genus ID and the ranking, just for neatness
        df = df.sort_values(["genus_id", "order_in_genus"])
        return df

    def _process_kraken(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Process the kraken classifications. Returns two dataframes:

        species_df: contains all the 'S' bacterial taxa and counts, with info added with respect to the parent genus. A
        confidence level (high or low) is applied based on filters.
        genus_df: contains all of 'G' bacterial taxa and counts, with added info about the species, such as the number
        that pass the filter.

        Also modifies the classifier_results dataframe by adding the parent taxonomy information.

        :params classifier_results: directly takes the kraken report as a pandas dataframe.
        :returns: species_df dataframe and genus_df dataframe. If there are no Kraken results (empty dataframe), then
        species and genus df are returned empty.
        """
        tp = self.taxaplease

        classifier_results = self.classifier_results.copy()  # Don't edit the original dataframe

        if classifier_results.empty:
            # If the classifier is empty, just return 2 empty dataframes for species and genus dfs, return early
            return pd.DataFrame(), pd.DataFrame()

        # Add parent taxonomy info
        classifier_results[["is_bacteria", "parent_record", "parent_human_readable"]] = classifier_results.apply(
            lambda row: self._get_parent_taxonomy(row["taxon_id"]),
            axis=1,
            result_type="expand",
        )

        # Get ALL genus level rows
        genus_df = classifier_results[
            (classifier_results["raw_rank"] == "G") & (classifier_results["is_bacteria"])
        ].copy()
        # Calculate proportion of reads at the species level compared to genus level
        genus_df["prop_species"] = 1 - (genus_df["count_direct"] / genus_df["count_descendants"])

        # get ALL species level rows for bacterial taxa
        species_df = classifier_results[
            (classifier_results["raw_rank"] == "S") & (classifier_results["is_bacteria"])
        ].copy()
        # add genus ID for the species
        species_df["genus_id"] = species_df["taxon_id"].apply(tp.get_genus_taxid)
        # count number of species
        total_species = species_df.groupby("genus_id").size().reset_index(name="total_species_identified").fillna(0)

        # Add the total species in a genus to the genus df
        genus_df = pd.merge(
            how="left", left=genus_df, right=total_species, left_on="taxon_id", right_on="genus_id"
        ).drop(columns=["genus_id"])

        # count number of species with >= 10 reads (READ_THRESHOLD value) and get the number of species that pass that
        # filter
        filtered_species = (
            species_df[(species_df["count_descendants"] >= self.thresholds["READ_THRESHOLD"])]
            .groupby("genus_id")
            .size()
            .reset_index(name="filtered_species_identified")
            .fillna(0)
        )
        genus_df = pd.merge(genus_df, filtered_species, left_on="taxon_id", right_on="genus_id", how="left").drop(
            columns=["genus_id"]
        )

        # merge genus info rows with species info - get the cols of interest and then merge
        genus_merge = genus_df[
            [
                "taxon_id",  # this is the same as genus_id
                "count_descendants",
                "total_species_identified",
                "filtered_species_identified",
            ]
        ].rename(columns={"count_descendants": "genus_level_reads"})

        species_df = (
            pd.merge(species_df, genus_merge, left_on="genus_id", right_on="taxon_id")
            .drop(columns=["taxon_id_y"])
            .rename(columns={"taxon_id_x": "taxon_id"})
        )

        # Add the percentage of reads for each species in the genus, plus to ranking/order of species in the genus
        species_df = self._add_pct_and_rank(species_df)

        # Add the confidence - high if pass, low if fail
        species_df["kraken_confidence"] = species_df.apply(
            lambda row: self._get_kraken_confidence_rating(
                count_descendants=row["count_descendants"],
                order_in_genus=row["order_in_genus"],
                pct_genus_reads=row["pct_genus_reads"],
            ),
            axis=1,
        )

        # Add the confidence to the genus_df - set to low for all of them:
        genus_df["kraken_confidence"] = "low"

        # Don't need the parent record dict in there...
        species_df = species_df.drop(columns=["parent_record"])
        genus_df = genus_df.drop(columns=["parent_record"])

        return species_df, genus_df

    def _get_kraken_results(self) -> tuple[str, dict, pd.DataFrame, pd.DataFrame]:
        """
        Get the headline result and the results from Kraken for bacteria.
        :param sample_id: string of climb_id.
        :param original_kraken_results: pandas dataframe containing the original kraken results from Scylla.
        :param kraken_thresholds_dict: dictionary containing the filter thresholds for kraken classifications.
        :param taxaplease_instance: instance of taxaplease, default is none and in this case will create a new instance.
        :return: tuple; headline_result (str), result (dict), kraken species result, kraken genus result.
        """

        kraken_species, kraken_genus = self._process_kraken()

        if kraken_species.empty:
            headline_result = "Kraken classified 0 bacterial species."
            results = {}
            return headline_result, results, kraken_species, kraken_genus

        # Add profiles to species
        kraken_species = profiler.add_profile_to_results(kraken_species, self.profiles_dict, self.taxaplease)

        # Add profiles to genus
        kraken_genus = profiler.add_profile_to_results(kraken_genus, self.profiles_dict, self.taxaplease)

        # Combine the species and the genera (in case there are genera without species calls)
        add_these_genera = kraken_genus[~kraken_genus["taxon_id"].isin(kraken_species["genus_id"])].copy()
        add_these_genera = add_these_genera.drop("prop_species", axis=1)
        combined_species_genera = pd.concat([kraken_species, add_these_genera], ignore_index=True)

        # Results is a dict of the profile taxa that were matched (genus and species), keep high and low confidence:
        unique_profile_taxa = combined_species_genera.sort_values(["profile", "kraken_confidence"])
        unique_profile_taxa = unique_profile_taxa[
            ["profile_taxon_match", "profile_taxon_id", "profile", "kraken_confidence"]
        ].drop_duplicates()

        results = unique_profile_taxa.dropna().reset_index(drop=True).to_dict(orient="index")

        #### Make headline result
        # Make headline result with format "ProfileA-high, ProfileA-low".
        profile_confidence = unique_profile_taxa.apply(
            lambda row: f"{row['profile']}-{row['kraken_confidence']}", axis=1
        )
        profile_confidence = sorted(set(profile_confidence.to_list()))
        total_kraken_species = kraken_species.shape[0]
        total_kraken_genera = add_these_genera.shape[0]

        headline_result = f"Kraken classified {total_kraken_species} bacterial species and {total_kraken_genera} genera; {', '.join(profile_confidence)}"

        return headline_result, results, kraken_species, kraken_genus

    # Make analysis tables
    def get_kraken_bacteria_analysis_table(self):
        """
        Pull together all the class attributes into the analysis table.
        """
        analysis_table, _ = create_analysis_fields(
            domain="bacteria",
            classifier="kraken",
            record_id=self.sample_id,
            thresholds=self.thresholds,
            headline_result=self.headline_result,
            results=self.results,
            server=self.server,
        )

        self.analysis_table = analysis_table

        return analysis_table

    def save_outputs_to_csv(self, results_dir: str | Path) -> None:
        """
        Save the final results to csv.
        :param filename: str, name of file to save to.
        :param results_dir: str or path to directory to save to.
        """
        species_unique_filename = f"{self.sample_id}_kraken_processed_species"
        write_df_to_csv(df=self.kraken_species_results, filename=species_unique_filename, results_dir=results_dir)

        genus_unique_filename = f"{self.sample_id}_kraken_processed_genera"
        write_df_to_csv(df=self.kraken_genus_results, filename=genus_unique_filename, results_dir=results_dir)


##########
# Sylph: #
##########
class SylphBacteria:
    """
    Class for parsing the bacterial classifications from Sylph.

    The main methods are:
    get_sylph_analysis_table - get the analysis table (returns instance of the onyx analysis helper class).
    save_outputs_to_csv - method to save the outputs (saved in instance attributes) to file (returns None)

    Atrributes:
    classifier_results - the original sylph outputs from Scylla.
    thresholds - dict; the thresholds to filter on.
    sample_id - str; climb-id
    profiles_dict - dict; lookup for profiles and their taxa.
    taxaplease - instance of taxaplease. Will instantiate a new taxaplease instance one if not given.
    server - str; database server.
    sylph_filtered_results - pd.DataFrame: All the taxa sylph identified for the sample, plus the confidence.
    headline_results - str; the main result, automatically generated to include the final number of taxa that were
     assigned high confidence.
    results - dict; the sylph_filtered_results dataframe filtered to high confidence species as a dict.
    analysis_table - oa.OnyxAnalysis; instance of the analysis table from the helper, containing all the relevant info.


    :param sample_id: str, climb-id
    :param original_classifier_df: pandas dataframe, the original results from scylla.
    :param sylph_bacteria_thresholds_dict: dict, containing the thresholds to filter for sylph.
    :param profiles_dict: dict containing profiles and their taxa.
    :param taxaplease_instance: instance of TaxaPlease class (optional)
    :param server: str, database server for Onyx to connect to - will be validated.
    """

    def __init__(
        self,
        sample_id: str,
        original_sylph_df: pd.DataFrame,
        sylph_bacteria_thresholds_dict: dict[str, int | float],
        profiles_dict: dict,
        taxaplease_instance: TaxaPlease | None,
        server: str = "server",
    ):
        """
        Create instance of SylphBacteria class, where arguments are attributes, and instance methods populate
        headline_result, result and sylph_processed_df.

        :param sample_id: str, climb-id
        :param original_classifier_df: pandas dataframe, the original results from scylla.
        :param sylph_bacteria_thresholds_dict: dict, containing the thresholds to filter for sylph.
        :param profiles_dict: dict containing profiles and their taxa.
        :param taxaplease_instance: instance of TaxaPlease class (optional).
        :param server: str, database server for Onyx to connect to - will be validated.
        """
        self.sylph: pd.DataFrame = original_sylph_df
        self.thresholds: dict = sylph_bacteria_thresholds_dict
        self.sample_id: str = sample_id
        self.profiles_dict: dict[str, dict[int, dict[str, str]]] = profiles_dict
        self.taxaplease: TaxaPlease = taxaplease_instance if taxaplease_instance else TaxaPlease()
        self.server: str = server

        self.headline_result: str
        self.result: dict
        self.sylph_processed_df: pd.DataFrame

        self.headline_result, self.results, self.sylph_processed_df = self._get_sylph_results()

    def _process_sylph_rank(self, row: pd.Series) -> tuple[int | None, str | None]:
        """
        Get the species taxon ID and the species name from row, using taxaplease.

        :params row: pd.series, row of a dataframe (used with an apply).
        :return: list of two; taxon ID and human-readable species name, or [None, None] if rank is not species or
        strain.
        """
        if row["taxon_id"] is None:
            logging.error("Taxon id is not found for %s" % (row))
            return None, None

        if row["taxon_rank"] == "species":
            return row["taxon_id"], row["human_readable"]
        elif row["taxon_rank"] == "strain":
            tp = self.taxaplease
            taxon_id = int(row["taxon_id"])
            species = tp.get_record(tp.get_species_taxid(taxon_id))  # type: ignore
            if species:
                return species["taxid"], species["name"]
            else:
                return None, None
        else:
            logging.info(
                "Taxon ID returned rank other than species or strain: %s, %s" % (row["taxon_id"], row["taxon_rank"])
            )
            return None, None

    def _get_sylph_confidence_rating(
        self,
        *,
        containment_index: float,
        effective_coverage: float,
    ) -> str:
        """
        Get the confidence rating for sylph results using thresholds.
        """
        if (
            containment_index >= self.thresholds["CONTAINMENT_INDEX_THRESHOLD"]
            and effective_coverage >= self.thresholds["EFFECTIVE_COVERAGE_THRESHOLD"]
        ):
            return "high"
        else:
            return "low"

    def _process_sylph(
        self,
    ) -> pd.DataFrame:
        """
        Process the sylph outputs and apply filters. Normalise to species level (sylph uses reference genomes which
        could be species or strain level, however the taxonomy parsing could lead to taxa at genus level), and add a
        confidence rating using the filters.
        :return: dataframe, with extra columns for filters. If there are no sylph results, return empty df.
        """

        tp = self.taxaplease

        # Get taxonomic rank of the sylph results
        sylph_df = self.sylph.copy()

        if sylph_df.empty:
            return pd.DataFrame()

        sylph_df["taxon_rank"] = sylph_df["taxon_id"].apply(
            lambda x: rec["rank"] if (rec := tp.get_record(x)) is not None else "Unknown"
        )

        sylph_df[["species_id", "species_human_readable"]] = sylph_df.apply(
            lambda x: self._process_sylph_rank(x), axis=1, result_type="expand"
        )

        sylph_df["genus_id"] = sylph_df["taxon_id"].apply(
            lambda x: tp.get_genus_taxid(x) if x is not None else None
        )  # Note that genus is not always the parent (unclassified, complexes)!

        sylph_df["cont_ind_eval"] = sylph_df["containment_index"].apply(eval)

        sylph_df["sylph_confidence"] = sylph_df.apply(
            lambda row: self._get_sylph_confidence_rating(
                containment_index=row["cont_ind_eval"], effective_coverage=row["effective_coverage"]
            ),
            axis=1,
        )
        return sylph_df

    def _get_sylph_results(self) -> tuple[str, dict, pd.DataFrame]:
        """
        Apply the filtering and add profiles, then get the headline result and the results from Sylph.

        Headline result: Sylph Classifier found taxa belonging to Profile 1A-low, profile1B-high...

        Result: a unique list of the profiles and profile taxa.

        :return: tuple; headline_result (str), result (dict), final table from sylph filtering (all the results with the
        added columns) to write to csv. If there are no sylph results, return headline result (automatically produced),
        results as empty dict and sylph_processed_df as empty dataframe.
        """
        sylph_processed_df = self._process_sylph()

        if sylph_processed_df.empty:
            headline_result = "Sylph classified 0 bacterial (or archaeal) taxa."
            results = {}
            return headline_result, results, sylph_processed_df

        # Add profiles
        sylph_processed_df = profiler.add_profile_to_results(sylph_processed_df, self.profiles_dict, self.taxaplease)

        # Make a table of the unique profiles with high confidence:
        unique_profile_taxa = sylph_processed_df.sort_values(["profile", "sylph_confidence"])
        unique_profile_taxa = unique_profile_taxa[
            ["profile_taxon_match", "profile_taxon_id", "profile", "sylph_confidence"]
        ].drop_duplicates()

        # Results is a dict of the profile taxa that were matched and only high confidence:
        results = (
            unique_profile_taxa.loc[unique_profile_taxa["sylph_confidence"] == "high"]
            .dropna()  # don't want to keep unknown profile taxa in the results table
            .reset_index(drop=True)
            .to_dict(orient="index")
        )
        # Make headline result with format "ProfileA (high, low)" with the profile if found, and high or low if found.

        #### Make headline result
        # Make headline result with format "ProfileA-high, ProfileA-low".
        profile_confidence = unique_profile_taxa.apply(
            lambda row: f"{row['profile']}-{row['sylph_confidence']}", axis=1
        )
        profile_confidence = sorted(set(profile_confidence.to_list()))

        total_sylph = sylph_processed_df.shape[0]
        headline_result = f"Sylph classified {total_sylph} bacterial (or archael) taxa; {', '.join(profile_confidence)}"

        return headline_result, results, sylph_processed_df

    def get_sylph_analysis_table(self):
        """
        Pull together all the class attributes into the analysis table.
        """
        analysis_table, _ = create_analysis_fields(
            domain="bacteria",
            classifier="sylph",
            record_id=self.sample_id,
            thresholds=self.thresholds,
            headline_result=self.headline_result,
            results=self.results,
            server=self.server,
        )

        self.analysis_table = analysis_table

        return analysis_table

    def save_outputs_to_csv(self, results_dir: str | Path) -> None:
        """
        Save the final results to csv.
        :param filename: str, name of file to save to.
        :param results_dir: str or path to directory to save to.
        """
        filename = f"{self.sample_id}_sylph_processed"
        write_df_to_csv(df=self.sylph_processed_df, filename=filename, results_dir=results_dir)
