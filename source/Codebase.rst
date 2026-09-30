########
Codebase
########

The codebase is Python based, and requires Profiler and TaxaPlease.
It also requires a configuation file which is provided but a user may define their own.
Furthermore, a clinical profile lookup json file is required to apply a profile to the taxa.
Finally a taxaplease database may be specified, however if none is supplied then TaxaPlease will download and create one.

The code base is split into bacteria and viruses, with general functionality available in utils and the commandline
application pulled together in main.

.. toctree::
   :maxdepth: 4

   bacteria
   viruses
   utils
