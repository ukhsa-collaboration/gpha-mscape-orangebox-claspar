import os

os.environ["CONFIG"] = "Placeholder config"
os.environ["ONYX_DOMAIN"] = "Placeholder domain"
os.environ["ONYX_TOKEN"] = "Placeholder token"


# def get_taxaplease_db(path_to_db):
#     tp = TaxaPlease()
#     tp = TaxaPlease(database=path_to_db)
#     tp.set_taxonomy_url("https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_2026-03-01.zip")


# @pytest.fixture(scope="session")
# def taxaplease_db(tmpdir_factory):

#     path_to_db = Path(tmpdir_factory.mktemp("data").join("taxa.db"))
#     get_taxaplease_db(path_to_db)
#     return path_to_db
