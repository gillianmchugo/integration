# interval files
## prepare and list interval files
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(readr)) install.packages("readr")
if (!require(stringr)) install.packages("stringr")
if (!require(tidyr)) install.packages("tidyr")

### list clumped files
clumped <- list.files("10_integration-analysis/01_input/04_intervals",
                      full.names = T,
                      pattern = "clumped")

### prepare associated interval files
for(i in clumped){
  readr::read_table(i) %>%
    dplyr::mutate(pos1 = BP - 1000000,
                  pos1 = dplyr::if_else(pos1 < 0,
                                        0,
                                        pos1),
                  pos2 = BP + 1000000) %>%
    dplyr::select(CHR,
                  pos1,
                  pos2) %>%
    readr::write_delim(paste0(stringr::str_remove(i,
                                                  "-plink-input.txt"),
                              "-intervals.txt"),
                       col_names = F)
}

### list interval files
interval_file <- list.files("10_integration-analysis/01_input/04_intervals",
                            pattern = ".clumped-intervals.txt")

interval_info <- data.frame(length = NA,
                            interval = NA)
for (i in interval_file){
  lines <- readr::read_delim(paste0("10_integration-analysis/01_input/04_intervals/",
                                    i),
                             col_names = F) %>%
    dplyr::count() %>%
    dplyr::pull()
  info <- data.frame(length = lines,
                     interval = stringr::str_remove(i,
                                                    pattern = ".clumped-intervals.txt"))
  interval_info <- dplyr::full_join(interval_info,
                                    info)
}
interval_info %>%
  tidyr::drop_na() %>%
  readr::write_csv("10_integration-analysis/01_input/04_intervals/interval-info.csv")