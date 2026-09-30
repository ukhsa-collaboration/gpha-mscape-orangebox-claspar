###########
🖱️ Usage 🖱️
###########

On the commandline, once installed, run:

.. code-block:: bash

    claspar --sample_id <ID> --output_dir <DIR> --server <server> --profiles_json <PATH-TO-CLINICAL_PROFILES.json>

Where `ID` is a valid sample ID and `DIR` is a valid output directory. This method includes an Onyx call, and therefore
requires all the necessary credentials and connectivity. There is a test profiles json in the tests/test_profile_tables directory.

However, it is possible to run ClasPar with a samplesheet (in 2x2 format, tab-seperated, with the columns 'climb-id' and
'full_Onyx_json'):

.. code-block:: bash

    claspar --sample_id <ID> --output_dir <DIR> --server server --samplesheet <TSV FILE> --profiles_json <PATH-TO-CLINICAL_PROFILES.json>

***NOTE**: sample ID and server is still required.*

*************
📃 Inputs: 📃
*************
+---------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Argument            | Required | Description                                                                                                                                                        |
+=====================+==========+====================================================================================================================================================================+
| --sample_id, -i     | Yes      | Climb-ID for sample                                                                                                                                                |
+---------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| --output_dir, -o    | Yes      | Path to directory where results will be saved to. The directory is made if it doesn't exist already.                                                               |
+---------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| --config, -c        | No       | Path to yaml file with filtering thresholds. Default can be seen in src/claspar/data/filter_thresholds.yaml                                                        |
+---------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| --server, -s        | Yes      | Must be one of: [server, synthscape]. Specify server code is being run on - helpful if developing on synthscape and running on server.                             |
+---------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| --samplesheet, -t   | No       | Path to samplesheet. Must be tsv, should have header 'full_Onyx_json' or 2 columns, 1 row.                                                                         |
+---------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| --profiles_json, -p | Yes      | Path to profiles table. Must be in .json format and have the required format. See Profiler https://github.com/ukhsa-collaboration/gpha-mscape-orangebox-profiler   |
+---------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| --log_file, -l      | No       | Path to log file. Default will be a file called '/sample-id/_claspar_/date-time/.log' in the output directory (where /date-time/ is the date and time of running). |
+---------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| --database_path, -d | No       | Path to database for TaxaPlease.                                                                                                                                   |
+---------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| --version, -v       | No       | print the version and exit.                                                                                                                                        |
+---------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+

***************
📤 Outputs: 📤
***************

There are a few outputs for a single sample, plus a log file. The log file appears in the same location as the outputs.

* ``<sample_id>_<data-time>_claspar.log`` - This is the log file and will contain helpful information and progress.

* ``<sample_id>_filtered_viral_aligner_results.csv`` - These are the viral aligner results filtered to those that pass the
  thresholds in the supplied config, plus the profiles for each taxa.

* ``<sample_id>_<classifier>_analysis_fields.json`` - This is the analysis table in json format, containing all the
  required analysis information. There is one per classifier.

* ``<sample_id>_kraken_processed_genera.csv`` - This is a table of all of the taxa at a genus level Kraken called. It has
  some additional columns with info about the species within the genus, and the profiles for these taxa at the genus
  level.

* ``<sample_id>__kraken_processed_species.csv`` - This is a table of all of the taxa at a species level Kraken called. It
  has some additional columns with info about the parent taxa (at genus level), total species in genus identified (and
  species that pass the filters), the proportion of total genus reads, the rank of that in its genus and the kraken
  confidence (high if it passes all of the filters, else low), and the profile for that taxa at the species level.

* ``<sample_id>_sylph_processed.csv``` - The sylph results after applying filtering and confidence (high if is passes all
  of the filters, else low), plus the profile for each taxa.

***********************
🔨 Troubleshooting: 🔨
***********************

* `claspar -h` returns list of help information for commands.
* The log should contain adequate information if there is an error or warning.
* There are unit tests available with test data. Run `pytest .` in the root of the repo after cloning. (You might need
  to have installed using 'Developer' settings.).
* The commandline entry point is in the `main()` function in main.py. There is a class per parser, split between
  bacteria.py and virus.py.
* package entry points are through `claspar.main.main`, which in turn calls `claspar.utils` (which houses some set up
  functionality like config parsing), `claspar.bacteria` (which houses the `KrakenBacteria` and `SylphBacteria` classes)
  and `claspar.virus` (which houses the 'VirusClasPar' class)
