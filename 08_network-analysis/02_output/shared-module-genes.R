# shared module genes
## extract genes shared between modules
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(readr)) install.packages("readr")
if (!require(stringr)) install.packages("stringr")

### read base network
base_network_sif <- readr::read_delim("08_network-analysis/01_input/base-network.sif",
                                      col_names = F) %>%
  dplyr::rename(gene1 = X1,
                interaction = X2,
                gene2 = X3)

### extract base network gene symbols
base_network_symbol <- dplyr::bind_rows(dplyr::select(base_network_sif,
                                                      gene1) %>%
                                          dplyr::rename(symbol = gene1),
                                        dplyr::select(base_network_sif,
                                                      gene2)%>%
                                          dplyr::rename(symbol = gene2)) %>%
  dplyr::distinct()

### set modules
modules <- c("MICRO_BL_34",
             "MICRO_LI_35",
             "MICRO_LN_35",
             "MICRO_SP_35",
             "RNA_LAGU_40",
             "RNA_BAOU_40",
             "RNA_NDAM_40",
             "RNA_BORG_40")

### extract gene symbols from modules
networks <- data.frame(symbol = NA, module = NA)
for (i in modules){
  sif <- readr::read_delim(paste0("08_network-analysis/02_output/",
                                  stringr::str_replace_all(i,
                                                           "_",
                                                           "-"),
                                  "-module-01.sif"),
                           col_names = F) %>%
    dplyr::rename(gene1 = X1,
                  interaction = X2,
                  gene2 = X3)
  symbols <- dplyr::bind_rows(dplyr::select(sif,
                                            gene1) %>%
                                dplyr::rename(symbol = gene1),
                              dplyr::select(sif,
                                            gene2)%>%
                                dplyr::rename(symbol = gene2)) %>%
    dplyr::distinct()
  test <- base_network_symbol %>%
    dplyr::mutate(module := dplyr::case_when(symbol %in% dplyr::pull(symbols) ~ stringr::str_replace_all(i,
                                                                                                         "_",
                                                                                                         " "),
                                             .default = NA))
  networks <- dplyr::full_join(networks,
                               test)
}

### identify and save shared gene symbols
shared_module_genes <- networks %>%
  tidyr::drop_na() %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(modules = toString(module),
                   number = dplyr::n()) %>%
  dplyr::filter(!stringr::str_detect(symbol,
                                     "BT.|IDB-")) %>%
  dplyr::filter(number >= 4) %>%
  dplyr::arrange(dplyr::desc(number)) %>%
  dplyr::select(symbol,
                modules) %>%
  readr::write_csv("08_network-analysis/02_output/shared-module-genes.csv")