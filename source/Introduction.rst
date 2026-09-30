############
Introduction
############

===============================================
👻 ClasPar 👻 - the friendly classifier parser
===============================================

ClasPar parses, filters and writes classifier results to analysis tables.
This includes taxonomy wrangling and clinical profile look ups for classified taxa.

This is a Python-based tool with a commandline interface and `claspar` package, set up in such a way as
to be run in the Orange-Box within a NextFlow workflow.

*It has been developed specifically for mSCAPE and the Orange Box.*

****************
What does it do?
****************

ClasPar parses, filters and writes results and analysis tables for classifiers. Specifically, it handles Kraken (bacteria),
Sylph (bacteria/archaea) and viral_aligner (viruses) outputs.

It applies filters using the config file in ``src/claspar/data/filter_thresholds.yaml`` by default (or a custom config
file can be defined).

It uses the OnyxAnalysisHelper to write all the results to json files, was well as all the results to csv files.

It uses the Profiler codebase (pinned to v2.0.0) to lookup classified taxa in the clinical profile tables.

It uses TaxaPlease to wrangle taxonomy (pinned to v2.2.1).
