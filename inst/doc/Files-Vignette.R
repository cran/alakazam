## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Set the file paths from inside the package directory
# These files are smaller versions of the example databases previously mentioned
changeo_file <- system.file("extdata", "example_changeo.tab.gz", package="alakazam")
airr_file <- system.file("extdata", "example_airr.tsv.gz", package="alakazam")

# Read in the data
db_changeo <- alakazam::readChangeoDb(changeo_file)
db_airr <- airr::read_rearrangement(airr_file)

## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
#  # Write the data to a tab-delimited file
#  alakazam::writeChangeoDb(db_changeo, "changeo.tsv")
#  airr::write_rearrangement(db_airr, "airr.tsv")

