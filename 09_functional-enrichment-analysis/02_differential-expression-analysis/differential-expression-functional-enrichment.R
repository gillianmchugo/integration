# differential expression functional enrichment
## functional enrichment of rna-seq differential expression results
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

### read input data
response_edger_13107 <- readr::read_delim("07_differential-expression-analysis/01_output/response_edger_13107_IN_IPA.txt")

### set contrast ids
rna_response_contrast_ids <- c("LAGU_20",
                               "LAGU_30",
                               "LAGU_40",
                               "NDAM_20",
                               "NDAM_30",
                               "NDAM_40",
                               "BAOU_20",
                               "BAOU_30",
                               "BAOU_40",
                               "BORG_20",
                               "BORG_30",
                               "BORG_40")

### extract all de genes
rna_de_genes <- data.frame(ID = NA)
for (i in 1:length(rna_response_contrast_ids)){
  rna_contrast_de_genes <- response_edger_13107 %>%
    dplyr::filter(!!rlang::sym(paste0("FDR_",
                                      rna_response_contrast_ids[i])) <=  0.05) %>%
    dplyr::distinct(ID,
                    .keep_all = T) %>%
    dplyr::select(ID)
  rna_de_genes <- dplyr::full_join(rna_de_genes,
                                   rna_contrast_de_genes)
}
rna_de_genes <- rna_de_genes %>%
  tidyr::drop_na() %>%
  readr::write_delim("09_functional-enrichment-analysis/02_differential-expression-analysis/rna-de-genes.csv")

### function to generate colour palette
smooth_rainbow <- khroma::colour("smooth rainbow")

### generate colour palette
pal <- smooth_rainbow(31,
                      range = c(0.1,
                                0.9))

### set background
background <- response_edger_13107 %>%
  dplyr::pull(ID)

### gprofiler functional analysis
go <- gprofiler2::gost(correction_method = "gSCS",
                       evcodes = T,
                       highlight = T,
                       organism = "btaurus",
                       sources = c("GO"),
                       query = dplyr::pull(rna_de_genes,
                                           ID))

### select top ten highlighted terms for each query
top_terms <- subset(go[["result"]],
                    highlighted == T) %>%
  dplyr::arrange(p_value) %>%
  dplyr::slice(1:10) %>%
  dplyr::mutate(top_term = T)

### add top term variable and save gprofiler results
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
                term_name_wrap = stringr::str_wrap(term_name, width = 10),
                query_label = factor(query,
                                     labels = paste0("'RNA'"))) %>%
  readr::write_csv("09_functional-enrichment-analysis/02_differential-expression-analysis/rna-gprofiler.csv")

### get numbers of terms
bp_terms <- go[["meta"]][["result_metadata"]][["GO:BP"]][["number_of_terms"]]
cc_terms <- go[["meta"]][["result_metadata"]][["GO:CC"]][["number_of_terms"]]
mf_terms <- go[["meta"]][["result_metadata"]][["GO:MF"]][["number_of_terms"]]
all_terms <- sum(bp_terms,
                 cc_terms,
                 mf_terms)

### plot gprofiler results
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
  ggplot2::continuous_scale(aesthetics = c("size", "point.size"),
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

### save gprofiler plot
png("09_functional-enrichment-analysis/02_differential-expression-analysis/rna-gprofiler.png",
    height = 6,
    res = 300,
    units = "in",
    width = 9)
print(plot)
dev.off()
pdf("09_functional-enrichment-analysis/02_differential-expression-analysis/rna-gprofiler.pdf",
    height = 6,
    width = 9)
print(plot)
dev.off()