# Claspar v2.1.1
Added onyx hash to methods.

## Changed
- bumped version of onyx analysis helper to 0.6.0.
- onyx hash is written into methods by the onyx analysis helper method.

Somewhat arbitrary version bump to try to fix docker container version issues.

---

# Claspar v2.1.0
Updated onyx query.

## Changed:
- Onyx query is now done once using the onyx analysis helper function. This returns the onyx versions which are fed into
the onyx analysis table object.
- tool versions are collated in main and fed into the analysis table creation function for each classifier.
- fixed the version of taxaplease
- bumped version of profiler (patch to fix taxaplease version ResolutionImpossible issue)
- docker publish dev as well as main.
- bumped version of onyx analysis helper to 0.6.0.

## Added:
- unit tests with patched queries.


---

# Claspar v2.0.1
2 small features added that modify the analysis table.
- added version of Profiler and the profiles table file used to the methods section of analysis tables
- changed the 'name' field in the analysis table to be more descriptive (now reads claspar-{classifier}-{domain} where classifier is viralaligner/sylph/kraken and domain is virus/bacteria).

---

# Claspar v2.0.0

BREAKING CHANGE - must now give the clinical profile tables file to claspar as argument.

- The clinical profile has been added to the claspar outputs. For every taxa identified by the three classifiers (including species and genera for Kraken), these are compared to the profile tables. Profile searching is hierarchical, which means if a species is identified by a classifier by only a taxon of higher taxonomic rank is present in a profile, this will be returned. Four columns are available for this -
'profile' - the profile the taxa belongs to.
'profile_taxon_match' - the name of the taxon the profile matches to.
'profile_taxon_id' - the taxon ID of the taxon matched in the profile.
'profile_rank' - the rank of the taxon matched in the profile.
- Analysis table outputs are created as json files, one per classifier. The headline result is a searchable string that contains the profiles and associated confidence identified in a sample. For example - "Viral Aligner classified 35 viral taxa; ProfileB-high, probileB-low, ProfileX-high, profileX-high".
- The 'results-metrics' in the analysis table is now a table of the profile taxa that have been seen, but not the specific taxa classified. (e.g. if E .coli, Salmonella and Klebsiella were seen in a sample, but only Enterobacteriaceae is in the profile table, just Enterobacteriaceae is reported in the table as one row.) The results-metrics will ONLY show High confidence for Sylph and High confidence for viral aligner, and high and low for kraken, with seperate rows for a taxon if both high and low confidence taxa were identified for that profile taxon.
- Kraken genera are included in the main result outputs. The 'results-metrics' table above will include any genera that Kraken identify at the genus level that have no species classified. These automatically get given the confidence level 'low' because they had no species identified. The genera csv however lists all taxa identified at the genus level, irrelevant of whether any species were found or not. It is possible to see whether species were found for a given genus using the column 'total_species_identified'.

Other changes:
- The mean read identity filter for viral aligner outputs (MEAN_READ_IDENTITY: 90) will be decreased from 90 to 80. This will inform whether a virus the viral aligner classified gets 'high' or 'low' confidence in the confidence column in the final output, and thus will be included in the analysis table 'results_metrics' (high) or not (low) - this is for the 'epi' outputs.
- Ranking species within a genus for Kraken taxa has been tweaked. There is a columm called 'order_in_genus' in the species csv, where 1 means 1st, 2 means 2nd etc. If two species have the same proportion in a genus, they both get the same place, with the next species getting the next place. So for a genus with 4 species with counts 15, 10, 10, 3 reads, they would get 1st, 2nd, 2nd, 3rd. All would pass the 'genus rank' (GENUS_RANK_THRESHOLD) filter (but not necessarily the 'proportion of genus' (GENUS_READ_PCT_THRESHOLD) filter).
- Claspar version is available on the command line using `--version` argument.
- Optional argument to supply a path to a TaxaPlease database.


Claspar v1.0.0
- parsed classifier outputs for kraken bacteria, sylph and viral aligner.
- analysis tables created.
