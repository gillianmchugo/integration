# integration input
## prepare inputs for integration analysis
### install required packages
if (!require(BiocManager)) install.packages("BiocManager")
if (!require(biomaRt)) BiocManager::install("biomaRt")
if (!require(dplyr)) install.packages("dplyr")
if (!require(readr)) install.packages("readr")
if (!require(rlang)) install.packages("rlang")
if (!require(tidyr)) install.packages("tidyr")

### prepare reference gene file
#### select ensembl data set
ensembl <- biomaRt::useEnsembl(biomart = "genes",
                               dataset="btaurus_gene_ensembl")

#### select ensembl attributes
biomaRt::getBM(attributes = c("chromosome_name",
                              "start_position",
                              "end_position",
                              "ensembl_gene_id",
                              "external_gene_name"),
               mart = ensembl) %>%
  readr::write_delim("10_integration-analysis/01_input/01_genes/genes.txt",
                     col_names = F)

### prepare reference snp files
readr::read_delim("01_data-sources/02_filtered-snp-data/hd.bim",
                  col_names = F) %>%
  dplyr::select(X1,
                X4) %>%
  readr::write_delim("10_integration-analysis/01_input/02_snps/hd-snps.txt",
                     col_names = F)

readr::read_delim("01_data-sources/02_filtered-snp-data/ld.bim",
                  col_names = F) %>%
  dplyr::select(X1,
                X4) %>%
  readr::write_delim("10_integration-analysis/01_input/02_snps/ld-snps.txt",
                     col_names = F)

readr::read_delim("03_filter-data/05_filter-snps/ld.bim",
                  col_names = F) %>%
  dplyr::select(X1,
                X4) %>%
  readr::write_delim("10_integration-analysis/01_input/02_snps/merged-ld-snps.txt",
                     col_names = F)

### prepare target gene set file
#### prepare end module genes
module_files <- list.files("08_network-analysis/02_output/")
end_module_files <- module_files[grepl("34|35|40",
                                       module_files)]
end_module_genes <- data.frame(ensembl = NA,
                               set = NA)
for (i in 1:length(end_module_files)){
  sif <- readr::read_delim(paste0("08_network-analysis/02_output/",
                                  end_module_files[i]),
                           col_names = F) %>%
    dplyr::rename(gene1 = X1,
                  interaction = X2,
                  gene2 = X3)
  symbol <- dplyr::bind_rows(dplyr::select(sif,
                                           gene1) %>%
                               dplyr::rename(symbol = gene1),
                             dplyr::select(sif,
                                           gene2)%>%
                               dplyr::rename(symbol = gene2)) %>%
    dplyr::distinct()
  ensembl <- gprofiler2::gconvert(query = symbol$symbol,
                                  organism = "btaurus") %>%
    dplyr::filter(stringr::str_detect(target_number,
                                      ".1")) %>%
    dplyr::distinct(name,
                    .keep_all = T) %>%
    dplyr::transmute(ensembl = target,
                     set = stringr::str_remove(end_module_files[i],
                                               ".sif"))
  end_module_genes <- dplyr::full_join(end_module_genes,
                                       ensembl)
}

#### write to targets file
end_module_genes %>%
  tidyr::drop_na() %>%
  readr::write_delim("10_integration-analysis/01_input/03_targets/end-module-targets.txt",
                     col_names = F)

### prepare associated interval files
#### function to prepare plink input
plink_input <- function(group, data = "hd", software = "elai", population = NA){
  column <- dplyr::case_match(population, NA ~ "",
                              .default = paste0(population, "_"))
  readr::read_csv(paste0("01_data-sources/09_local-ancestry-results/",
                         software,
                         "-",
                         data,
                         "-",
                         stringr::str_replace_all(group,
                                                  "_",
                                                  "-"),
                         "-hybrids.csv")) %>%
    dplyr::mutate(p_a = pnorm(get(paste0(column,
                                         "z_a")),
                              lower.tail = F),
                  p_b = pnorm(get(paste0(column,
                                         "z_b")),
                              lower.tail = F),
                  p_c = pnorm(get(paste0(column,
                                         "z_c")),
                              lower.tail = F)) %>%
    readr::write_delim(paste0("10_integration-analysis/01_input/04_intervals/",
                              software,
                              "-",
                              data,
                              "-",
                              stringr::str_replace_all(column,
                                                       "_",
                                                       "-"),
                              stringr::str_replace_all(group,
                                                       "_",
                                                       "-"),
                              "-hybrids-plink-input.txt"))
}

#### function to prepare plink input for samples with microarray data
plink_input_microarray <- function(population, data = "hd", software = "elai"){
  readr::read_csv(paste0("06_local-ancestry-analysis/02_output/07_previous-results/",
                         software,
                         "-",
                         data,
                         "-",
                         population,
                         "-microarray.csv")) %>%
    dplyr::mutate(p_a = pnorm(z_a,
                              lower.tail = F),
                  p_b = pnorm(z_b,
                              lower.tail = F),
                  p_c = pnorm(z_c,
                              lower.tail = F)) %>%
    readr::write_delim(paste0("10_integration-analysis/01_input/04_intervals/",
                              software,
                              "-",
                              data,
                              "-",
                              population,
                              "-microarray-plink-input.txt"))
}

#### function to prepare plink input for samples with expression data
plink_input_expression <- function(population, expression = T){
  readr::read_csv(paste0("06_local-ancestry-analysis/02_output/",
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
                         ".csv")) %>%
    dplyr::mutate(p_a = pnorm(z_a,
                              lower.tail = F),
                  p_b = pnorm(z_b,
                              lower.tail = F),
                  p_c = pnorm(z_c,
                              lower.tail = F)) %>%
    readr::write_delim(paste0("10_integration-analysis/01_input/04_intervals/",
                              population,
                              "-expression",
                              dplyr::case_when(expression == T ~ "",
                                               expression == F ~ "-all"),
                              "-plink-input.txt"))
}

#### set groups
groups <- c("european",
            "selected_european",
            "trypanotolerant_african",
            "selected_trypanotolerant_african",
            "trypanosusceptible_african")

#### generate plink input files from previous results for groups
lapply(groups,
       plink_input)
lapply(groups,
       plink_input,
       data = "ld")

lapply(groups,
       plink_input,
       software = "mosaic")

#### generate plink input files for populations with expression samples
plink_input("trypanotolerant_african",
            population = "NDAM")
plink_input("trypanosusceptible_african",
            population = "BORA")

plink_input("trypanotolerant_african",
            data = "ld",
            population = "NDAM")
plink_input("trypanosusceptible_african",
            data = "ld",
            population = "BORA")

plink_input("trypanotolerant_african",
            population = "NDAM",
            software = "mosaic")
plink_input("trypanosusceptible_african",
            population = "BORA",
            software = "mosaic")

plink_input_microarray("NDAM")
plink_input_microarray("BORA")

plink_input_microarray("NDAM",
                       data = "ld")
plink_input_microarray("BORA",
                       data = "ld")

plink_input_expression("BAOU")
plink_input_expression("BORA")
plink_input_expression("BORG")
plink_input_expression("FULA")
plink_input_expression("LAGU")
plink_input_expression("NDAM")

plink_input_expression("BAOU",
                       expression = F)
plink_input_expression("BORA",
                       expression = F)
plink_input_expression("BORG",
                       expression = F)
plink_input_expression("FULA",
                       expression = F)
plink_input_expression("LAGU",
                       expression = F)
plink_input_expression("NDAM",
                       expression = F)