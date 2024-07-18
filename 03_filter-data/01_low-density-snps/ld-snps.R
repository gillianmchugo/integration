# ld snps
## write list of ld snps
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(readr)) install.packages("readr")

### read ld data and write list of ld snps
readr::read_table("03_filter-data/01_low-density-snps/berthier-widde.map",
                  col_names = F) %>%
  dplyr::select(X2) %>%
  readr::write_delim("03_filter-data/01_low-density-snps/ld-snps.txt",
                     col_names = F)