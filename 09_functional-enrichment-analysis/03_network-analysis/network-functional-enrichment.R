# network functional enrichment
## functional enrichment of network results
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggh4x)) install.packages("ggh4x")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggrepel)) install.packages("ggrepel")
if (!require(ggtext)) install.packages("ggtext")
if (!require(gprofiler2)) install.packages("gprofiler2")
if (!require(khroma)) install.packages("khroma")
if (!require(readr)) install.packages("readr")
if (!require(scales)) install.packages("scales")
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

### convert base network gene symbols to ensembl ids and save
base_network_ensembl <- gprofiler2::gconvert(query = base_network_symbol$symbol,
                                             organism = "btaurus") %>%
  dplyr::filter(stringr::str_detect(target_number,
                                    ".1")) %>%
  dplyr::distinct(name,
                  .keep_all = T) %>%
  dplyr::select(target) %>%
  readr::write_csv("09_functional-enrichment-analysis/03_network-analysis/01_input/base-network.csv")

### read base network ensembl ids
base_network_ensembl <- readr::read_csv("09_functional-enrichment-analysis/03_network-analysis/01_input/base-network.csv")

### function for gprofiler functional analysis of network results
network_gprofiler <- function(data = "microarray", module_number = 1){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palettes
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  pal2 <- khroma::colour("sunset")(11)
  # set background
  background <- base_network_ensembl %>%
    dplyr::pull(target)
  modules <- get(paste0(data,
                        "_modules"))
  input <- list()
  for (i in modules){
    sif <- readr::read_delim(paste0("08_network-analysis/02_output/",
                                    i,
                                    "-module-0",
                                    module_number,
                                    ".sif"),
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
    # set input
    input[[i]] <- gprofiler2::gconvert(query = symbol$symbol,
                                       organism = "btaurus") %>%
      dplyr::filter(stringr::str_detect(target_number,
                                        ".1")) %>%
      dplyr::distinct(name,
                      .keep_all = T) %>%
      dplyr::pull(target)
    # save input data
    input[[i]] %>%
      as.data.frame() %>%
      readr::write_csv(paste0("09_functional-enrichment-analysis/03_network-analysis/01_input/",
                              i,
                              "-module-0",
                              module_number,
                              ".csv"),
                       col_names = F)
  }
  # gprofiler functional analysis
  go <- gprofiler2::gost(correction_method = "gSCS",
                         custom_bg = background,
                         domain_scope = "custom_annotated",
                         evcodes = T,
                         highlight = T,
                         organism = "btaurus",
                         sources = c("GO"),
                         query = input)
  # select top ten highlighted terms for each query
  top_terms <- subset(go[["result"]], highlighted == T) %>%
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
                  query_label = paste0("'",
                                       stringr::str_replace_all(query,
                                                                "-",
                                                                " "),
                                       "'"),
                  rna_query_label = factor(query_label,
                                           levels = c("'RNA LAGU 40'",
                                                      "'RNA BAOU 40'",
                                                      "'RNA NDAM 40'",
                                                      "'RNA BORG 40'")))  %>%
    readr::write_csv(paste0("09_functional-enrichment-analysis/03_network-analysis/02_output/",
                            data,
                            "-modules-0",
                            module_number,
                            "-gprofiler.csv"))
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
    ggh4x::facet_grid2(dplyr::case_when(data == "microarray" ~ ggplot2::vars(query_label),
                                        data == "rna" ~ ggplot2::vars(rna_query_label)),
                       ggplot2::vars(source_label),
                       drop = F,
                       labeller = ggplot2::label_parsed,
                       scales = "free_x",
                       strip = ggh4x::strip_vanilla(clip = "off"),
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
    ggplot2::scale_colour_manual(values = c("RNA-BAOU-40" = pal2[[1]],
                                            "RNA-BORG-40" = pal[4],
                                            "RNA-LAGU-40" = pal[8],
                                            "RNA-NDAM-40" = pal[5],
                                            "MICRO-BL-34" = pal[28],
                                            "MICRO-LI-35" = pal[17],
                                            "MICRO-LN-35" = pal[9],
                                            "MICRO-SP-35" = pal[21])) +
    ggplot2::scale_fill_manual(values = ggplot2::alpha(c("RNA-BAOU-40" = pal2[[1]],
                                                         "RNA-BORG-40" = pal[4],
                                                         "RNA-LAGU-40" = pal[8],
                                                         "RNA-NDAM-40" = pal[5],
                                                         "MICRO-BL-34" = pal[28],
                                                         "MICRO-LI-35" = pal[17],
                                                         "MICRO-LN-35" = pal[9],
                                                         "MICRO-SP-35" = pal[21]),
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
                   axis.title.y = ggtext::element_markdown(),
                   legend.margin = ggplot2::margin(0,
                                                   0,
                                                   0,
                                                   0),
                   legend.position = "top",
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggtext::element_markdown(colour = "black",
                                                         size = ggplot2::rel(1)))
  # save gprofiler plot
  png(paste0("09_functional-enrichment-analysis/03_network-analysis/03_figures/",
             data,
             "-modules-0",
             module_number,
             "-gprofiler.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("09_functional-enrichment-analysis/03_network-analysis/03_figures/",
             data,
             "-modules-0",
             module_number,
             "-gprofiler.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

### set microarray modules
microarray_modules <- c("MICRO-BL-34",
                        "MICRO-LI-35",
                        "MICRO-LN-35",
                        "MICRO-SP-35")

### set rna-seq modules
rna_modules <- c("RNA-BAOU-40",
                 "RNA-BORG-40",
                 "RNA-LAGU-40",
                 "RNA-NDAM-40")

### apply function
network_gprofiler()
network_gprofiler(data = "rna")

### function to generate colour palette
smooth_rainbow <- khroma::colour("smooth rainbow")

### generate colour palette
pal <- smooth_rainbow(31,
                      range = c(0.1,
                                0.9))

### micro li 35 module
#### read and set input data
input <- readr::read_csv("09_functional-enrichment-analysis/03_network-analysis/01_input/MICRO-LI-35-module-01.csv",
                         col_names = F)

#### gprofiler functional analysis
go <- gprofiler2::gost(correction_method = "gSCS",
                       domain_scope = "annotated",
                       evcodes = T,
                       highlight = T,
                       organism = "btaurus",
                       sources = c("GO"),
                       query = dplyr::pull(input,
                                           X1))

#### select top ten highlighted terms for each query
top_terms <- subset(go[["result"]],
                    highlighted == T) %>%
  dplyr::arrange(p_value) %>%
  dplyr::slice(1:10) %>%
  dplyr::mutate(top_term = T)

#### add top term variable and save gprofiler results
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
                                     labels = paste0("'MICRO LI 35'"))) %>%
  readr::write_csv("09_functional-enrichment-analysis/03_network-analysis/02_output/MICRO-LI-35-module-01-gprofiler-no-background.csv")

#### get numbers of terms
bp_terms <- go[["meta"]][["result_metadata"]][["GO:BP"]][["number_of_terms"]]
cc_terms <- go[["meta"]][["result_metadata"]][["GO:CC"]][["number_of_terms"]]
mf_terms <- go[["meta"]][["result_metadata"]][["GO:MF"]][["number_of_terms"]]
all_terms <- sum(bp_terms, cc_terms, mf_terms)

#### plot gprofiler results
plot <- ggplot2::ggplot(go[["result"]],
                        ggplot2::aes(x = source_order,
                                     y = plot_p_value)) +
  ggplot2::geom_point(ggplot2::aes(colour = source,
                                   fill = source,
                                   size = intersection_ratio),
                      shape = 21,
                      stroke = 0.5) +
  ggh4x::facet_grid2(ggplot2::vars(query_label),
                     ggplot2::vars(source_label),
                     drop = F,
                     labeller = ggplot2::label_parsed,
                     scales = "free_x",
                     strip = ggh4x::strip_vanilla(clip = "off"),
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
  ggplot2::scale_colour_manual(values = c("GO:BP" = pal[26],
                                          "GO:CC" = pal[17],
                                          "GO:MF" = pal[28])) +
  ggplot2::scale_fill_manual(values = ggplot2::alpha(c("GO:BP" = pal[26],
                                                       "GO:CC" = pal[17],
                                                       "GO:MF" = pal[28]),
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
                 axis.title.y = ggtext::element_markdown(),
                 legend.margin = ggplot2::margin(0,
                                                 0,
                                                 0,
                                                 0),
                 legend.position = "top",
                 strip.background = ggplot2::element_blank(),
                 strip.text = ggtext::element_markdown(colour = "black",
                                                       size = ggplot2::rel(1)))

#### save gprofiler plot
png("09_functional-enrichment-analysis/03_network-analysis/03_figures/MICRO-LI-35-module-01-gprofiler-no-background.png",
    height = 6,
    res = 300,
    units = "in",
    width = 9)
print(plot)
dev.off()
pdf("09_functional-enrichment-analysis/03_network-analysis/03_figures/MICRO-LI-35-module-01-gprofiler-no-background.pdf",
    height = 6,
    width = 9)
print(plot)
dev.off()

### base network
#### set input data
input <- base_network_ensembl

#### gprofiler functional analysis
go <- gprofiler2::gost(correction_method = "gSCS",
                       domain_scope = "annotated",
                       evcodes = T,
                       highlight = T,
                       organism = "btaurus",
                       sources = c("GO"),
                       query = dplyr::pull(input,
                                           target))

#### select top ten highlighted terms for each query
top_terms <- subset(go[["result"]],
                    highlighted == T) %>%
  dplyr::arrange(p_value) %>%
  dplyr::slice(1:10) %>%
  dplyr::mutate(top_term = T)

#### add top term variable and save gprofiler results
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
                                     labels = paste0("'BASE'"))) %>%
  readr::write_csv("09_functional-enrichment-analysis/03_network-analysis/02_output/base-network-gprofiler.csv")

#### get numbers of terms
bp_terms <- go[["meta"]][["result_metadata"]][["GO:BP"]][["number_of_terms"]]
cc_terms <- go[["meta"]][["result_metadata"]][["GO:CC"]][["number_of_terms"]]
mf_terms <- go[["meta"]][["result_metadata"]][["GO:MF"]][["number_of_terms"]]
all_terms <- sum(bp_terms, cc_terms, mf_terms)

#### plot gprofiler results
plot <- ggplot2::ggplot(go[["result"]],
                        ggplot2::aes(x = source_order,
                                     y = plot_p_value)) +
  ggplot2::geom_point(ggplot2::aes(colour = source,
                                   fill = source,
                                   size = intersection_ratio),
                      shape = 21,
                      stroke = 0.5) +
  ggh4x::facet_grid2(ggplot2::vars(query_label),
                     ggplot2::vars(source_label),
                     drop = F,
                     labeller = ggplot2::label_parsed,
                     scales = "free_x",
                     strip = ggh4x::strip_vanilla(clip = "off"),
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
                           force = 100,
                           force_pull = 0,
                           lineheight = 0.7,
                           max.iter = 1000000000,
                           max.overlaps = Inf,
                           max.time = 1000,
                           min.segment.length = 0,
                           segment.size = 0.25,
                           size = ggplot2::rel(3)) +
  ggplot2::labs(size = "Intersection ratio",
                y = "-Log<sub>10</sub>*P*<sub>adj.</sub>") +
  ggplot2::guides(colour = "none",
                  fill = "none",
                  point.size = "none",
                  size = guide_legend(nrow = 1,
                                      order = 1)) +
  ggplot2::scale_colour_manual(values = c("GO:BP" = pal[26],
                                          "GO:CC" = pal[17],
                                          "GO:MF" = pal[28])) +
  ggplot2::scale_fill_manual(values = ggplot2::alpha(c("GO:BP" = pal[26],
                                                       "GO:CC" = pal[17],
                                                       "GO:MF" = pal[28]),
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
  theme(axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggtext::element_markdown(),
        legend.margin = ggplot2::margin(0,
                                        0,
                                        0,
                                        0),
        legend.position = "top",
        strip.background = ggplot2::element_blank(),
        strip.text = ggtext::element_markdown(colour = "black",
                                              size = ggplot2::rel(1)))

#### save gprofiler plot
png("09_functional-enrichment-analysis/03_network-analysis/03_figures/base-network-gprofiler.png",
    height = 6,
    res = 300,
    units = "in",
    width = 9)
print(plot)
dev.off()
pdf("09_functional-enrichment-analysis/03_network-analysis/03_figures/base-network-gprofiler.pdf",
    height = 6,
    width = 9)
print(plot)
dev.off()