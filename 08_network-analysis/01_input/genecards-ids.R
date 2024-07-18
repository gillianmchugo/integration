# genecards ids
## prepare genecards ids
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(gprofiler2)) install.packages("gprofiler2")
if (!require(readr)) install.packages("readr")

### read and filter genecards results
genecards <- readr::read_csv("01_data-sources/04_genecards/GeneCards-SearchResults.csv") %>%
  dplyr::filter(`Relevance score` >= 1.75)

### convert gene symbols to ensembl ids
genecards_ensembl <- gprofiler2::gconvert(query = genecards$`Gene Symbol`,
                                          organism = "btaurus") %>%
  dplyr::filter(stringr::str_detect(target_number,
                                    ".1")) %>%
  dplyr::distinct(name,
                  .keep_all = T) %>%
  dplyr::select(target) %>%
  readr::write_delim("04_network-analysis/01_input/genecards-ensembl.txt",
                     col_names = F,
                     delim = "\t")