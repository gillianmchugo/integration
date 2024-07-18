# local ancestry functional enrichment
## functional enrichment of local ancestry results
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggh4x)) install.packages("ggh4x")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggrepel)) install.packages("ggrepel")
if (!require(ggtext)) install.packages("ggtext")
if (!require(gprofiler2)) install.packages("gprofiler2")
if (!require(khroma)) install.packages("khroma")
if (!require(parallel)) install.packages("parallel")
if (!require(readr)) install.packages("readr")
if (!require(scales)) install.packages("scales")
if (!require(stringr)) install.packages("stringr")

### function to get genes for each chromosome
chromosome_genes <- function(chromosome, window = 1){
  # set snp positions
  positions <- readr::read_csv("06_local-ancestry-analysis/02_output/01_BAOU/elai-ld-BAOU.csv") %>%
    dplyr::filter(chr == chromosome) %>%
    dplyr::mutate(pos1 = pos - (window*1000000),
                  pos1 = if_else(pos1 < 0,
                                 0,
                                 pos1),
                  pos2 = pos + (window*1000000),
                  input = paste0(chr,
                                 ":",
                                 pos1,
                                 ":",
                                 pos2))
  # get snp positions
  snps <- positions %>%
    dplyr::pull(input)
  # convert to genes
  genes <- gprofiler2::gconvert(query = snps,
                                organism = "btaurus")
  # merge genes with snps and save file
  positions %>%
    dplyr::select(rs,
                  pos,
                  chr,
                  pos1,
                  pos2,
                  input) %>%
    dplyr::full_join(genes) %>%
    readr::write_csv(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/01_input/01_genes/elai-ld-",
                            window,
                            "mb-genes-",
                            stringr::str_pad(chromosome,
                                             2,
                                             pad = "0",
                                             side = "left"),
                            ".csv"))
}

### function to combine chromosome genes into genome genes
genome_genes <- function(window = 1){
  # get chromosome 1
  join <- readr::read_csv(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/01_input/01_genes/elai-ld-",
                                 window,
                                 "mb-genes-01.csv"))
  # add chromosomes 2 - 29
  for(i in 2:29){
    join <- join %>%
      dplyr::full_join(readr::read_csv(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/01_input/01_genes/elai-ld-",
                                              window,
                                              "mb-genes-",
                                              stringr::str_pad(i,
                                                               2,
                                                               pad = "0",
                                                               side = "left"),
                                              ".csv")))
  }
  # save file
  join %>%
    readr::write_csv(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/01_input/01_genes/elai-ld-",
                            window,
                            "mb-genes.csv"))
}

### function to generate elai functional analysis input file
functional_input <- function(population, window = 1, expression = F){
  # read input data
  input <- readr::read_csv(paste0("06_local-ancestry-analysis/02_output/",
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
  # merge analysis with genes and save to file
  background <- readr::read_csv(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/01_input/01_genes/elai-ld-",
                                       window,
                                       "mb-genes.csv")) %>%
    dplyr::select(c(rs,
                    pos,
                    chr,
                    target,
                    name)) %>%
    dplyr::full_join(input) %>%
    dplyr::select(c(rs,
                    pos,
                    chr,
                    target,
                    name,
                    mean_a,
                    mean_b,
                    mean_c,
                    z_a,
                    z_b,
                    z_c)) %>%
    dplyr::group_by(target,
                    name) %>%
    dplyr::summarize(mean_a = mean(mean_a),
                     mean_b = mean(mean_b),
                     mean_c = mean(mean_c),
                     z_a = mean(z_a),
                     z_b = mean(z_b),
                     z_c = mean(z_c)) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(p_value_a = pnorm(z_a,
                                    lower.tail = F),
                  p_value_b = pnorm(z_b,
                                    lower.tail = F),
                  p_value_c = pnorm(z_c,
                                    lower.tail = F),
                  name = dplyr::na_if(name, target)) %>%
    dplyr::relocate(name) %>%
    dplyr::rename(ensembl = target,
                  symbol = name) %>%
    readr::write_csv(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/01_input/02_analyses/elai-ld-",
                            population,
                            dplyr::case_when(expression == T ~ "-expression-",
                                             expression == F ~ "-"),
                            window,
                            "mb.csv"))
}

### function for gprofiler functional analysis of local ancestry results
gprofiler <- function(population, window = 1, zscore = 2, expression = F){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palette
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # read input data
  input <- readr::read_csv(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/01_input/02_analyses/elai-ld-",
                                  population,
                                  dplyr::case_when(expression == T ~ "-expression-",
                                                   expression == F ~ "-"),
                                  window,
                                  "mb.csv"))
  # set background
  background <- input %>%
    dplyr::pull(ensembl)
  # set query genes
  a <- input %>%
    dplyr::filter(z_a >= zscore) %>%
    dplyr::select(ensembl)
  b <- input %>%
    dplyr::filter(z_b >= zscore) %>%
    dplyr::select(ensembl)
  c <- input %>%
    dplyr::filter(z_c >= zscore) %>%
    dplyr::select(ensembl)
  # gprofiler functional analysis
  go <- gprofiler2::gost(correction_method = "gSCS",
                         custom_bg = background,
                         domain_scope = "custom_annotated",
                         evcodes = T,
                         highlight = T,
                         organism = "btaurus",
                         sources = c("GO"),
                         query = list("a" = a,
                                      "b" = b,
                                      "c" = c))
  # select top ten highlighted terms for each query
  top_terms <- subset(go[["result"]],
                      highlighted == T) %>%
    dplyr::group_by(query) %>%
    dplyr::arrange(p_value) %>%
    dplyr::slice(1:10) %>%
    dplyr::mutate(top_term = T) %>%
    dplyr::ungroup()
  # add top term variable and save gprofiler results
  go[["result"]] <- dplyr::full_join(go[["result"]],
                                     top_terms) %>%
    dplyr::mutate(intersection_ratio = intersection_size/term_size,
                  plot_p_value = dplyr::case_when(-log10(p_value) > 16 ~ 16,
                                                  .default = -log10(p_value)),
                  source_label = factor(source,
                                        labels = c("GO:BP" = "'GO:Biological Process'",
                                                   "GO:CC" = "'GO:Cellular<br>Component'",
                                                   "GO:MF" = "'GO:Molecular Function'"),
                                        levels = c("GO:BP",
                                                   "GO:CC",
                                                   "GO:MF")),
                  term_name_wrap = stringr::str_wrap(term_name,
                                                     width = 10),
                  query_label = factor(query,
                                       labels = c("a" = "'European'~italic('Bos taurus')",
                                                  "b" = "'African'~italic('Bos taurus')",
                                                  "c" = "italic('Bos indicus')"),
                                       levels = c("a",
                                                  "b",
                                                  "c"))) %>%
    readr::write_csv(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/02_output/gprofiler-elai-ld-",
                            population,
                            dplyr::case_when(expression == T ~ "-expression",
                                             expression == F ~ ""),
                            ".csv"))
  # get numbers of terms
  bp_terms <- go[["meta"]][["result_metadata"]][["GO:BP"]][["number_of_terms"]]
  cc_terms <- go[["meta"]][["result_metadata"]][["GO:CC"]][["number_of_terms"]]
  mf_terms <- go[["meta"]][["result_metadata"]][["GO:MF"]][["number_of_terms"]]
  all_terms <- sum(bp_terms,
                   cc_terms,
                   mf_terms)
  # plot gprofiler results
  plot <- ggplot2::ggplot(go[["result"]],
                          ggplot2::aes(x = source_order,
                                       y = plot_p_value)) +
    ggplot2::geom_point(ggplot2::aes(colour = query,
                                     fill = query,
                                     size = intersection_ratio),
                        shape = 21,
                        stroke = 0.5) +
    ggh4x::facet_grid2(ggplot2::vars(query_label),
                       ggplot2::vars(source_label),
                       drop = F,
                       labeller = ggplot2::label_parsed,
                       scales = "free_x",
                       strip = strip_vanilla(clip = "off"),
                       switch = "x") +
    ggplot2::geom_point(data = subset(go[["result"]],
                                      top_term == T),
                        ggplot2::aes(size = intersection_ratio),
                        shape = 21,
                        stroke = 0.5) +
    ggplot2::continuous_scale(aesthetics = c("size",
                                             "point.size"),
                              breaks = c(0.25,
                                         0.5,
                                         0.75,
                                         1),
                              labels = c(0.25,
                                         0.5,
                                         0.75,
                                         1),
                              limits = c(0,
                                         1),
                              scale_name = "size",
                              palette = scales::area_pal()) +
    ggrepel::geom_text_repel(data = subset(go[["result"]],
                                           top_term == T),
                             ggplot2::aes(label = term_name_wrap,
                                          point.size = intersection_ratio),
                             bg.color = "white",
                             bg.r = 0.1,
                             box.padding = 0.5,
                             force = 20,
                             force_pull = 0,
                             lineheight = 0.7,
                             max.iter = 1000000000,
                             max.overlaps = Inf,
                             max.time = 10,
                             min.segment.length = 0,
                             segment.size = 0.25,
                             size = ggplot2::rel(3)) +
    ggplot2::labs(size = "Intersection ratio",
                  y = "-Log<sub>10</sub>*P*<sub>adj.</sub>") +
    ggplot2::guides(colour = "none",
                    fill = "none",
                    point.size = "none",
                    size = ggplot2::guide_legend(nrow = 1,
                                                 order = 1)) +
    ggplot2::scale_colour_manual(values = c("a" = pal[17],
                                            "b" = pal[9],
                                            "c" = pal[28])) +
    ggplot2::scale_fill_manual(values = ggplot2::alpha(c("a" = pal[17],
                                                         "b" = pal[9],
                                                         "c" = pal[28]),
                                                       0.8)) +
    ggh4x::scale_x_facet(COL == 1,
                         breaks = c(0,
                                    cc_terms,
                                    cc_terms*2,
                                    cc_terms*3,
                                    cc_terms*4,
                                    cc_terms*5,
                                    cc_terms*6,
                                    cc_terms*7),
                         minor_breaks = NULL,
                         limits = c(0,
                                    bp_terms)) +
    ggh4x::scale_x_facet(COL == 2,
                         minor_breaks = NULL,
                         breaks = c(0,
                                    cc_terms),
                         limits = c(0,
                                    cc_terms)) +
    ggh4x::scale_x_facet(COL == 3,
                         minor_breaks = NULL,
                         breaks = c(0,
                                    cc_terms,
                                    cc_terms*2),
                         limits = c(0,
                                    mf_terms)) +
    ggplot2::scale_y_continuous(breaks = c(0,
                                           4,
                                           8,
                                           12,
                                           16),
                                expand = c(0,
                                           0),
                                labels = c(0,
                                           4,
                                           8,
                                           12,
                                           expression("">=16)),
                                limits = c(0,
                                           17)) +
    ggh4x::force_panelsizes(cols = c(bp_terms/all_terms,
                                     cc_terms/all_terms,
                                     mf_terms/all_terms),
                            respect = F) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::ggtext::element_markdown(),
                   legend.margin = ggplot2::margin(0,
                                                   0,
                                                   0,
                                                   0),
                   legend.position = "top",
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggtext::element_markdown(colour = "black",
                                                         size = ggplot2::rel(1)))
  # save gprofiler plot
  png(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/03_figures/gprofiler-elai-ld-",
             population,
             dplyr::case_when(expression == T ~ "-expression",
                              expression == F ~ ""),
             ".png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/03_figures/gprofiler-elai-ld-",
             population,
             dplyr::case_when(expression == T ~ "-expression",
                              expression == F ~ ""),
             ".pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

### apply functions to data
#### get genes
parallel::mclapply(1:29,
                   chromosome_genes,
                   mc.cores = 15)
genome_genes()

#### set list of analyses
populations <- c("BAOU",
                 "BORA",
                 "BORG",
                 "FULA",
                 "LAGU",
                 "NDAM")

#### functional analysis input
parallel::mclapply(populations,
                   functional_input,
                   mc.cores = 15)

parallel::mclapply(populations,
                   functional_input,
                   expression = T,
                   mc.cores = 15)

#### gprofiler functional analysis
gprofiler("BAOU")
gprofiler("BORA")
gprofiler("BORG")
gprofiler("FULA")
gprofiler("LAGU")
gprofiler("NDAM")

gprofiler("BAOU",
          expression = T)
gprofiler("BORA",
          expression = T)
gprofiler("BORG",
          expression = T)
gprofiler("FULA",
          expression = T)
gprofiler("LAGU",
          expression = T)
gprofiler("NDAM",
          expression = T)