# elai results
## extract elai results
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(parallel)) install.packages("parallel")
if (!require(readr)) install.packages("readr")
if (!require(scales)) install.packages("scales")
if (!require(stringr)) install.packages("stringr")
if (!require(tibble)) install.packages("tibble")

### function to read chromosome data and extract mean ancestry results for samples with expression data
elai_chromosome_data_expression <- function(chromosome, population, expression = T){
  # read mean local ancestry dosage file
  ps21 <- readr::read_table(paste0("06_local-ancestry-analysis/02_output/",
                                   dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                                                    population == "BORA" ~ "02_BORA",
                                                    population == "BORG" ~ "03_BORG",
                                                    population == "FULA" ~ "04_FULA",
                                                    population == "LAGU" ~ "05_LAGU",
                                                    population == "NDAM" ~ "06_NDAM"),
                                   "/output/ld-",
                                   population,
                                   "-",
                                   chromosome,
                                   ".ps21.txt"),
                            col_names = F) %>%
    dplyr::mutate(across(everything(),
                         ~./2))
  # read snp information file
  snpinfo <- readr::read_table(paste0("06_local-ancestry-analysis/02_output/",
                                      dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                                                       population == "BORA" ~ "02_BORA",
                                                       population == "BORG" ~ "03_BORG",
                                                       population == "FULA" ~ "04_FULA",
                                                       population == "LAGU" ~ "05_LAGU",
                                                       population == "NDAM" ~ "06_NDAM"),
                                      "/output/ld-",
                                      population,
                                      "-",
                                      chromosome,
                                      ".snpinfo.txt")) %>%
    dplyr::mutate(chr = chromosome)
  # set order
  order <- c(1,
             3,
             2)
  # extract ancestries
  a <- t(ps21[seq(order[1],
                  length(ps21),
                  3)])
  b <- t(ps21[seq(order[2],
                  length(ps21),
                  3)])
  c <- t(ps21[seq(order[3],
                  length(ps21),
                  3)])
  # extract sample ids
  fam <- readr::read_delim("03_filter-data/05_filter-snps/ld.fam",
                           col_names = F,
                           show_col_types = F) %>%
    dplyr::rename(population_id = X1,
                  id = X2) %>%
    dplyr::filter(population_id == population)
  # add sample ids to column names
  colnames(a) <- fam$id
  colnames(b) <- fam$id
  colnames(c) <- fam$id
  # filter sample ids
  expression_ids <- fam %>%
    dplyr::mutate(expression = stringr::str_detect(id,
                                                   dplyr::case_when(expression == T ~ "berthier|BORA_mchugo|NDAM_mchugo",
                                                                    expression == F ~ "_"))) %>%
    dplyr::filter(expression == T) %>%
    dplyr::pull(id)
  # extract samples with expression data
  a <- a[, expression_ids]
  b <- b[, expression_ids]
  c <- c[, expression_ids]
  # calculate mean ancestry
  mean_a <- apply(a,
                  1,
                  mean)
  mean_b <- apply(b,
                  1,
                  mean)
  mean_c <- apply(c,
                  1,
                  mean)
  # ensure correct length
  mean_a <- mean_a[1:nrow(snpinfo)]
  mean_b <- mean_b[1:nrow(snpinfo)]
  mean_c <- mean_c[1:nrow(snpinfo)]
  # combine ancestries and save file
  tibble::tibble(snpinfo,
                 mean_a,
                 mean_b,
                 mean_c) %>%
    readr::write_csv(paste0("06_local-ancestry-analysis/02_output/",
                            dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                                             population == "BORA" ~ "02_BORA",
                                             population == "BORG" ~ "03_BORG",
                                             population == "FULA" ~ "04_FULA",
                                             population == "LAGU" ~ "05_LAGU",
                                             population == "NDAM" ~ "06_NDAM"),
                            "/elai-ld-",
                            population,
                            "-",
                            stringr::str_pad(chromosome,
                                             2,
                                             pad = "0",
                                             side = "left"),
                            dplyr::case_when(expression == T ~ "-expression",
                                             expression == F ~ ""),
                            ".csv"))
}

### function to combine chromosomes into genome data for samples with expression data
elai_genome_data_expression <- function(population, expression = T){
  # get chromosome 1
  join <- readr::read_csv(paste0("06_local-ancestry-analysis/02_output/",
                                 dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                                                  population == "BORA" ~ "02_BORA",
                                                  population == "BORG" ~ "03_BORG",
                                                  population == "FULA" ~ "04_FULA",
                                                  population == "LAGU" ~ "05_LAGU",
                                                  population == "NDAM" ~ "06_NDAM"),
                                 "/elai-ld-",
                                 population,
                                 "-01",
                                 dplyr::case_when(expression == T ~ "-expression",
                                                  expression == F ~ ""),
                                 ".csv"))
  # add chromosomes 2 - 29
  for(i in 2:29){
    join <- join %>%
      dplyr::full_join(readr::read_csv(paste0("06_local-ancestry-analysis/02_output/",
                                              dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                                                               population == "BORA" ~ "02_BORA",
                                                               population == "BORG" ~ "03_BORG",
                                                               population == "FULA" ~ "04_FULA",
                                                               population == "LAGU" ~ "05_LAGU",
                                                               population == "NDAM" ~ "06_NDAM"),
                                              "/elai-ld-",
                                              population,
                                              "-",
                                              stringr::str_pad(i,
                                                               2,
                                                               pad = "0",
                                                               side = "left"),
                                              dplyr::case_when(expression == T ~ "-expression",
                                                               expression == F ~ ""),
                                              ".csv")))
  }
  # format snp positions, add row numbers and z-scores and save file
  join %>%
    dplyr::mutate(row_num = dplyr::row_number(),
                  z_a = as.numeric(scale(mean_a)),
                  z_b = as.numeric(scale(mean_b)),
                  z_c = as.numeric(scale(mean_c))) %>%
    readr::write_csv(paste0("06_local-ancestry-analysis/02_output/",
                            dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                                             population == "BORA" ~ "02_BORA",
                                             population == "BORG" ~ "03_BORG",
                                             population == "FULA" ~ "04_FULA",
                                             population == "LAGU" ~ "05_LAGU",
                                             population == "NDAM" ~ "06_NDAM"),
                            "/elai-ld-",
                            population,
                            dplyr::case_when(expression == T ~ "-expression",
                                             expression == F ~ ""),
                            ".csv"))
}

### apply functions to data
#### read elai chromosome data and extract mean ancestry results
##### samples with expression data
parallel::mclapply(1:29,
                   elai_chromosome_data_expression,
                   population = "BAOU",
                   mc.cores = 15)
parallel::mclapply(1:29,
                   elai_chromosome_data_expression,
                   population = "BORA",
                   mc.cores = 15)
parallel::mclapply(1:29,
                   elai_chromosome_data_expression,
                   population = "BORG",
                   mc.cores = 15)
parallel::mclapply(1:29,
                   elai_chromosome_data_expression,
                   population = "FULA",
                   mc.cores = 15)
parallel::mclapply(1:29,
                   elai_chromosome_data_expression,
                   population = "LAGU",
                   mc.cores = 15)
parallel::mclapply(1:29,
                   elai_chromosome_data_expression,
                   population = "NDAM",
                   mc.cores = 15)

##### all samples
parallel::mclapply(1:29,
                   elai_chromosome_data_expression,
                   population = "BAOU",
                   expression = F,
                   mc.cores = 15)
parallel::mclapply(1:29,
                   elai_chromosome_data_expression,
                   population = "BORA",
                   expression = F,
                   mc.cores = 15)
parallel::mclapply(1:29,
                   elai_chromosome_data_expression,
                   population = "BORG",
                   expression = F,
                   mc.cores = 15)
parallel::mclapply(1:29,
                   elai_chromosome_data_expression,
                   population = "FULA",
                   expression = F,
                   mc.cores = 15)
parallel::mclapply(1:29,
                   elai_chromosome_data_expression,
                   population = "LAGU",
                   expression = F,
                   mc.cores = 15)
parallel::mclapply(1:29,
                   elai_chromosome_data_expression,
                   population = "NDAM",
                   expression = F,
                   mc.cores = 15)

#### combine elai chromosome data into genome data for all hybrids
##### samples with expression data
elai_genome_data_expression("BAOU")
elai_genome_data_expression("BORA")
elai_genome_data_expression("BORG")
elai_genome_data_expression("FULA")
elai_genome_data_expression("LAGU")
elai_genome_data_expression("NDAM")

##### all samples
elai_genome_data_expression("BAOU",
                            expression = F)
elai_genome_data_expression("BORA",
                            expression = F)
elai_genome_data_expression("BORG",
                            expression = F)
elai_genome_data_expression("FULA",
                            expression = F)
elai_genome_data_expression("LAGU",
                            expression = F)
elai_genome_data_expression("NDAM",
                            expression = F)