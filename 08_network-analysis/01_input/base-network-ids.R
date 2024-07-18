# base network ids
## prepare base network
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(readr)) install.packages("readr")
if (!require(stringr)) install.packages("stringr")

### read, prepare and save base network
readr::read_delim("01_data-sources/05_innatedb/temp_batch_05385395257358292.sif",
                  col_names = F) %>%
  dplyr::transmute(gene1 = stringr::str_split_i(X1,
                                                " ",
                                                1),
                   interaction = X2,
                   gene2 = stringr::str_split_i(X3,
                                                " ",
                                                1)) %>% 
  readr::write_delim("04_network-analysis/01_input/base-network.sif",
                     col_names = F,
                     delim = "\t")