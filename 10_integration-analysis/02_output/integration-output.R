# integration output
## format output of integration analysis
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(readr)) install.packages("readr")
if (!require(stringr)) install.packages("stringr")

### list output files
files <- list.files("10_integration-analysis/02_output/01_input/04_intervals",
                    full.names = T)

### write list of output files
list.files("10_integration-analysis/02_output/01_input/04_intervals") %>%
  as.data.frame() %>%
  readr::write_csv("10_integration-analysis/02_output/analyses.csv")

### check for significant results
results <- c()
for(i in files){
  test <- readr::read_file(i) %>%
    stringr::str_detect("no significant results at p < 0.05")
  results <- c(results,
               dplyr::case_when(test == T ~ NA,
                                test == F ~ i))
}

### copy files with significant results to output folder
file.copy(from = results,
          to = "10_integration-analysis/02_output/02_results")