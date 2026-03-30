#!/usr/bin/env python3

"""
ClasPar - parse classifiers, apply filters and create analysis tables. Also run profiler to create profiles for sample.
"""

import argparse
import logging
import sys
from datetime import datetime
from importlib import resources
from pathlib import Path

from profiler import get_profiles as profiler
from taxaplease import TaxaPlease

from claspar import __version__, bacteria, virus
from claspar.utils import ClasParError, get_input_data, read_config_file, read_samplesheet, setup_outdir

today = datetime.today().strftime("%Y-%m-%d")


# Arg parse setup
def get_args():
    parser = argparse.ArgumentParser(
        prog="claspar",
        description=f"""
        ClasPar: the friendly classifier parser that parses, filters and writes classifier results to analysis tables.
        Version = {__version__}
        """,
    )
    parser.add_argument("--sample_id", "-i", dest="sample_id", type=str, required=True, help="Climb-ID for sample.")
    parser.add_argument(
        "--output_dir",
        "-o",
        dest="output_dir",
        type=str,
        required=True,
        help="Path to directory where results will be saved to. Directory will be created if it does not exist.",
    )
    parser.add_argument(
        "--config",
        "-c",
        dest="config",
        type=str,
        required=False,
        help="Optional - Path to yaml file with filtering thresholds",
    )
    parser.add_argument(
        "--server",
        "-s",
        dest="server",
        type=str,
        required=True,
        choices=["mscape", "synthscape"],
        help="Specify server code is being run on - helpful if developing on synthscape and running on server",
    )
    parser.add_argument(
        "--profiles-table",
        "-p",
        dest="profile_table_spreadsheet_path",
        type=Path,
        required=True,
        help="Path to profile tables spreadsheet. Must be in xlsx, with tabs named as the profiles to be assigned.",
    )
    parser.add_argument(
        "--database_path", "-d", type=str, required=False, help="Optional - Path to a database file for TaxaPlease."
    )
    parser.add_argument(
        "--samplesheet",
        "-t",
        dest="samplesheet_path",
        type=str,
        required=False,
        help="Optional - Path to samplesheet. Must be tsv, should have header 'full_Onyx_json' or 2 columns, 1 row.",
    )
    parser.add_argument(
        "--log-file",
        "-l",
        dest="log_file",
        type=str,
        required=False,
        help=(
            """
            Optional - Path to log file. Default will be a file called '/sample-id/_claspar_/date-time/.log' in the
            output directory (where sample_id is the climb-id for the sample and /date-time/ is the date and time of
            running).
            """
        ),
    )
    parser.add_argument("--version", "-v", action="version", version=f"%(prog)s - version {__version__}")

    return parser


# Logger set up
def set_up_logger(log_filepath):
    """
    Set up logger, set to append mode so logs from older runs are not overwritten.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s")

    out_handler = logging.FileHandler(log_filepath, mode="a")
    out_handler.setFormatter(formatter)
    logger.addHandler(out_handler)

    return logger


# Main function
def main():
    """Run ClasPar - parse the classifiers (class instance of each classifier), and run profiler
    to add clinical profiles to all taxa. Write out to files.
    Must exit cleanly on design - for integration with OB.
    """

    try:
        #########
        # Setup #
        #########

        # set up arg parser:
        parser = get_args()

        # If no args provided, print help and exit
        if len(sys.argv) == 1:
            parser.print_help()
            sys.exit(0)

        # Get the args
        args = parser.parse_args()

        # Set up output dir:
        args.output_dir = Path(args.output_dir)
        setup_outdir(args.output_dir)

        # Set up log file:
        log_file = (
            Path(args.log_file)
            if args.log_file
            else Path(args.output_dir) / f"{args.sample_id}_{today}_claspar_log.txt"
        )
        set_up_logger(log_file)

        # Set up thresholds:
        threshold_dict = {}
        # Use default filtering thresholds yaml file if custom file is not supplied
        if not args.config:
            with resources.as_file(resources.files("claspar.data").joinpath("filter_thresholds.yaml")) as config_file:
                threshold_dict, exit_codes = read_config_file(config_file)
                if any(exit_codes):
                    # logging happens in the function
                    # --> Exit if any issues with the config:
                    return 1

            logging.info(
                "No custom filtering thresholds yaml file specified, using default parameters from included file."
            )

        else:
            logging.info("Reading the filtering thresholds from custom yaml file provided: %s" % (args.config))

            # Read in filtering thresholds from yaml file
            try:
                threshold_dict, config_exit_codes = read_config_file(args.config)
                # If any of the exit_codes are 1, then some filters were missing - need to exit.
                # These are logged in the check_filters function:
                # --> Exit if any issues with the config:
                if any(config_exit_codes):
                    # logging happens in the function
                    return 1
            # --> Exit if config file not found:
            except FileNotFoundError:
                logging.error("Specified filtering thresholds yaml file %s not found, exiting program." % (args.config))
                return 1

        # Set up data here - samplesheet or onyx:
        if args.samplesheet_path:
            dataframes, samplesheet_exitcode = read_samplesheet(args.samplesheet_path)
            if samplesheet_exitcode == 1:
                # --> Exit if issues with samplesheet parsing:
                logging.error("Cannot parse samplesheet provided, exiting.")
                return 1

        else:
            # Set up data needed (query Onyx once here)
            dataframes, onyx_exitcode = get_input_data(args.sample_id, args.server)
            if onyx_exitcode == 1:
                # --> Exit if issues with Onyx:
                logging.error("Exiting due to issues with Onyx.")
                return 1

        # Unpack dataframes into variables:
        try:
            viral_aligner_input_df, sylph_input_df, classifier_calls_df = dataframes
        except ValueError as e:
            # --> Exit if cannot unpack dataframes:
            logging.error("Cannot unpack the dataframes into the variables: %s" % (e))  # noqa
            return 1

        # Set up profile lookup dict:
        profiles_dict = profiler.make_profiles_dict(args.profile_table_spreadsheet_path)

        tp = TaxaPlease(database=args.database_path) if args.database_path else TaxaPlease()

        ####################
        # The Actual Thing #
        ####################

        # Bacteria -

        ##########
        # Kraken #
        ##########

        # Add in rest of code including logging messages:
        logging.info("Parsing the Kraken Bacteria classifications...")

        # Instantiate the kraken bacterial parser class
        kraken_bacteria_parser = bacteria.KrakenBacteria(
            sample_id=args.sample_id,
            original_classifier_df=classifier_calls_df,
            kraken_bacteria_thresholds_dict=threshold_dict["kraken_bacterial_filters"],
            profiles_dict=profiles_dict,
            taxaplease_instance=tp,
            server=args.server,
        )

        # Get the analysis table:
        kraken_bacterial_analysis_table = kraken_bacteria_parser.get_kraken_bacteria_analysis_table()

        # All good so far, let's write analysis table to json:
        kraken_bacteria_json_path = (
            Path(args.output_dir) / f"{args.sample_id}.claspar-krakenbacteria.analysis_fields.json"
        )
        kraken_bacterial_analysis_table.write_analysis_to_json(result_file=kraken_bacteria_json_path)
        logging.info(
            "Kraken Bacterial Classifications for Onyx analysis fields written to file %s" % (kraken_bacteria_json_path)
        )

        # Write files to csv:
        kraken_bacteria_parser.save_outputs_to_csv(args.output_dir)
        logging.info("All Processed Kraken bacterial species and genera data written to csv in %s" % (args.output_dir))
        logging.info("Finished parsing Kraken Bacterial results.")

        #########
        # Sylph #
        #########

        # Add in rest of code including logging messages:
        logging.info("Parsing the Sylph classifications...")

        # Instantiate the kraken bacterial parser class
        sylph_parser = bacteria.SylphBacteria(
            sample_id=args.sample_id,
            original_sylph_df=sylph_input_df,
            sylph_bacteria_thresholds_dict=threshold_dict["sylph_filters"],
            profiles_dict=profiles_dict,
            taxaplease_instance=tp,
            server=args.server,
        )

        # Get the analysis table:
        sylph_analysis_table = sylph_parser.get_sylph_analysis_table()

        # All good so far, let's write analysis table to json:
        sylph_json_path = Path(args.output_dir) / f"{args.sample_id}.claspar-sylph.analysis_fields.json"
        sylph_analysis_table.write_analysis_to_json(result_file=sylph_json_path)

        logging.info("Sylph Bacterial Classifications for Onyx analysis fields written to file %s" % (sylph_json_path))

        # Write files to csv:
        sylph_parser.save_outputs_to_csv(args.output_dir)
        logging.info("All processed Sylph data written to csv in %s" % (args.output_dir))

        logging.info("Finished parsing Sylph results.")

        ###########
        # Viruses #
        ###########

        # Add in rest of code including logging messages:
        logging.info("Parsing the viral aligner classifications.")

        viral_aligner = virus.VirusClasPar(
            sample_id=args.sample_id,
            original_viral_aligner_df=viral_aligner_input_df,
            virus_thresholds_dict=threshold_dict["viral_aligner_filters"],
            profiles_dict=profiles_dict,
            taxaplease_instance=tp,
            server=args.server,
        )

        viral_aligner_analysis_table = viral_aligner.get_virus_analysis_table()

        # All good so far, let's write analysis table to json:
        viral_aligner_json_path = Path(args.output_dir) / f"{args.sample_id}.claspar-viralaligner.analysis_fields.json"

        viral_aligner_analysis_table.write_analysis_to_json(result_file=viral_aligner_json_path)

        logging.info("Viral Aligner Onyx analysis fields written to file %s", viral_aligner_json_path)

        # Save csv files to go to s3:
        viral_aligner.save_outputs_to_csv(results_dir=args.output_dir)
        logging.info("Viral Aligner filtered data written to csv in %s" % (args.output_dir))
        logging.info("Finished parsing Viral Aligner results.")

        ########
        # End! #
        ########

    except ClasParError:
        return 1

    except Exception as e:
        msg = "Unhandled error: %s." % (e)
        logging.exception(msg, exc_info=True)
        return 1

    return 0


# Run
if __name__ == "__main__":
    sys.exit(main())
