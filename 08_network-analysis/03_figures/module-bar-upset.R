# module bar upset
## module bar and upset plots
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(khroma)) install.packages("khroma")
if (!require(patchwork)) install.packages("patchwork")
if (!require(readr)) install.packages("readr")
if (!require(rlang)) install.packages("rlang")
if (!require(stringr)) install.packages("stringr")
if (!require(tidyr)) install.packages("tidyr")

### bar plot
#### read base network
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
modules <- c("RNA_BAOU_40",
             "RNA_BORG_40",
             "RNA_LAGU_40",
             "RNA_NDAM_40",
             "MICRO_BL_34",
             "MICRO_LI_35",
             "MICRO_LN_35",
             "MICRO_SP_35")

### extract gene symbols from modules
networks <- data.frame(symbol = NA)
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
    dplyr::mutate(!!rlang::sym(paste0(i)) := dplyr::case_when(symbol %in% dplyr::pull(symbols) ~ T,
                                                              .default = F))
  networks <- dplyr::full_join(networks,
                               test)
}

#### get data for bar plot
bar <- data.frame(module = NA, number = NA, category = NA)
for (i in modules){
  module_shared <- data.frame(module = stringr::str_replace_all(i,
                                                                "_",
                                                                " "),
                              total = networks %>%
                                tidyr::drop_na() %>%
                                dplyr::filter(!!rlang::sym(i) == T) %>%
                                dplyr::count() %>%
                                dplyr::pull(),
                              unique = networks %>%
                                tidyr::drop_na() %>%
                                dplyr::filter(!!rlang::sym(i) == T & if_all(modules[! modules %in% i], ~ . %in% F)) %>%
                                dplyr::count() %>%
                                dplyr::pull(),
                              category = "Shared") %>%
    dplyr::mutate(number = total - unique) %>%
    dplyr::select(module,
                  number,
                  category)
  module_unique <- data.frame(module = stringr::str_replace_all(i,
                                                                "_",
                                                                " "),
                              number =  networks %>%
                                tidyr::drop_na() %>%
                                dplyr::filter(!!rlang::sym(i) == T & if_all(modules[! modules %in% i], ~ . %in% F)) %>%
                                dplyr::count() %>%
                                dplyr::pull(),
                              category = stringr::str_replace_all(i,
                                                                  "_",
                                                                  " "))
  bar <- dplyr::full_join(bar,
                          module_shared)
  bar <- dplyr::full_join(bar,
                          module_unique)
}

#### function to generate colour palette
smooth_rainbow <- khroma::colour("smooth rainbow")

#### generate colour palettes
pal <- smooth_rainbow(31,
                      range = c(0.1,
                                0.9))

pal2 <- khroma::colour("sunset")(11)

#### draw bar plot
plot <- bar %>%
  tidyr::drop_na() %>%
  ggplot2::ggplot(ggplot2::aes(x = factor(module,
                                          levels = c("MICRO BL 34",
                                                     "MICRO LI 35",
                                                     "MICRO LN 35",
                                                     "MICRO SP 35",
                                                     "RNA LAGU 40",
                                                     "RNA BAOU 40",
                                                     "RNA NDAM 40",
                                                     "RNA BORG 40")),
                               y = number,
                               fill = factor(category,
                                             levels = c("MICRO BL 34",
                                                        "MICRO LI 35",
                                                        "MICRO LN 35",
                                                        "MICRO SP 35",
                                                        "RNA LAGU 40",
                                                        "RNA BAOU 40",
                                                        "RNA NDAM 40",
                                                        "RNA BORG 40",
                                                        "Shared")))) +
  ggplot2::geom_col() +
  ggplot2::labs(x = "Module",
                y = "Number of Genes") +
  ggplot2::scale_fill_manual(values = c(pal[28],
                                        pal[17],
                                        pal[9],
                                        pal[21],
                                        pal[8],
                                        pal2[[1]],
                                        pal[5],
                                        pal[4],
                                        "#666666")) +
  ggplot2::scale_y_continuous(expand = c(0,
                                         0),
                              limits = c(0,
                                         2500)) +
  ggplot2::theme_light() +
  ggplot2::theme(legend.title = ggplot2::element_blank(),
                 legend.position = "top")

#### save bar plot
png("08_network-analysis/03_figures/03_bar/module-bar.png",
    height = 6,
    res = 300,
    units = "in",
    width = 9)
print(plot)
dev.off()
pdf("08_network-analysis/03_figures/03_bar/module-bar.pdf",
    height = 6,
    width = 9)
print(plot)
dev.off()

### upset plot
#### get data for upset plot
upset <- networks %>%
  tidyr::drop_na() %>%
  dplyr::select(!symbol)

#### function to generate colour palette
smooth_rainbow <- khroma::colour("smooth rainbow")

#### generate colour palettes
pal <- smooth_rainbow(31,
                      range = c(0.1,
                                0.9))

pal2 <- khroma::colour("sunset")(11)

#### set labels
upset_labels <- c("RNA_BAOU_40" = "RNA BAOU 40",
                  "RNA_BORG_40" = "RNA BORG 40",
                  "RNA_LAGU_40" = "RNA LAGU 40",
                  "RNA_NDAM_40" = "RNA NDAM 40",
                  "MICRO_BL_34" = "MICRO BL 34",
                  "MICRO_LI_35" = "MICRO LI 35",
                  "MICRO_LN_35" = "MICRO LN 35",
                  "MICRO_SP_35" = "MICRO SP 35")

#### set queries
upset_queries <- list(ComplexUpset::upset_query(set = "MICRO_BL_34",
                                                fill = pal[28],
                                                only_components = "overall_sizes"),
                      ComplexUpset::upset_query(intersect = "MICRO_BL_34",
                                                fill = pal[28],
                                                only_components = "Intersection size"),
                      ComplexUpset::upset_query(set = "MICRO_LI_35",
                                                fill = pal[17],
                                                only_components = "overall_sizes"),
                      ComplexUpset::upset_query(intersect = "MICRO_LI_35",
                                                fill = pal[17],
                                                only_components = "Intersection size"),
                      ComplexUpset::upset_query(set = "MICRO_LN_35",
                                                fill = pal[9],
                                                only_components = "overall_sizes"),
                      ComplexUpset::upset_query(intersect = "MICRO_LN_35",
                                                fill = pal[9],
                                                only_components = "Intersection size"),
                      ComplexUpset::upset_query(set = "MICRO_SP_35",
                                                fill = pal[21],
                                                only_components = "overall_sizes"),
                      ComplexUpset::upset_query(intersect = "MICRO_SP_35",
                                                fill = pal[21],
                                                only_components = "Intersection size"),
                      ComplexUpset::upset_query(set = "RNA_BORG_40",
                                                fill = pal[4],
                                                only_components = "overall_sizes"),
                      ComplexUpset::upset_query(intersect = "RNA_BORG_40",
                                                fill = pal[4],
                                                only_components = "Intersection size"),
                      ComplexUpset::upset_query(set = "RNA_NDAM_40",
                                                fill = pal[5],
                                                only_components = "overall_sizes"),
                      ComplexUpset::upset_query(intersect = "RNA_NDAM_40",
                                                fill = pal[5],
                                                only_components = "Intersection size"),
                      ComplexUpset::upset_query(set = "RNA_BAOU_40",
                                                fill = pal2[[1]],
                                                only_components = "overall_sizes"),
                      ComplexUpset::upset_query(set = "RNA_LAGU_40",
                                                fill = pal[8],
                                                only_components = "overall_sizes"),
                      ComplexUpset::upset_query(intersect = "RNA_LAGU_40",
                                                fill = pal[8],
                                                only_components = "Intersection size"))

#### draw upset plot
upset_plot <- ComplexUpset::upset(upset,
                                  intersect = colnames(upset),
                                  encode_sets = F,
                                  guide = "over",
                                  min_degree = 1,
                                  n_intersections = 20,
                                  name = "Intersections",
                                  width_ratio = 0.25,
                                  labeller = ggplot2::as_labeller(upset_labels),
                                  queries = upset_queries,
                                  base_annotations = list('Intersection size' = ComplexUpset::intersection_size(mapping = ggplot2::aes(fill = "#000000"),
                                                                                                                linewidth = 1,
                                                                                                                width = 0.8,
                                                                                                                text = list(size = ggplot2::rel(3))) +
                                                            ggplot2::scale_y_continuous(expand = c(0,
                                                                                                   0),
                                                                                        limits = c(0,
                                                                                                   950)) +
                                                            ggplot2::scale_fill_identity() +
                                                            ggplot2::ylab("Number of Genes\n in Common") ),
                                  matrix = ComplexUpset::intersection_matrix(geom = ggplot2::geom_point(size = 1.5)),
                                  set_sizes = ComplexUpset::upset_set_size(geom = ggplot2::geom_bar(linewidth = 1,
                                                                                                    width = 0.75)) +
                                    scale_y_reverse(expand = c(0,
                                                               0),
                                                    limits = c(2500,
                                                               0)) +
                                    ylab("Total Number of Genes\nin Module"),
                                  stripes = ComplexUpset::upset_stripes(geom = ggplot2::geom_segment(alpha = 0.5,
                                                                                                     linewidth = 6.8),
                                                                        colors = c(pal[28],
                                                                                   pal[17],
                                                                                   pal[9],
                                                                                   pal[8],
                                                                                   pal2[[1]],
                                                                                   pal[4],
                                                                                   pal[5],
                                                                                   pal[21])),
                                  themes = list('Intersection size' = ggplot2::theme_light() +
                                                  ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                                 axis.ticks.x = ggplot2::element_blank(),
                                                                 axis.title.x = ggplot2::element_blank()),
                                                intersections_matrix = ggplot2::theme_light() +
                                                  ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                                 axis.ticks = ggplot2::element_blank(),
                                                                 axis.title.y = ggplot2::element_blank(),
                                                                 axis.text.y = ggplot2::element_text(family = "mono",
                                                                                                     margin = ggplot2::margin(l = -40)),
                                                                 legend.position = "none"),
                                                overall_sizes = ggplot2::theme_light() +
                                                  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                                                 axis.ticks.y = ggplot2::element_blank(),
                                                                 axis.title.y = ggplot2::element_blank(),
                                                                 plot.margin = ggplot2::margin(r = 13))))

#### get upset data for legend plot
upset_data <- upset_labels %>%
  as.data.frame() %>%
  dplyr::rename(label = ".") %>%
  dplyr::mutate(label = factor(label,
                               levels = c("MICRO BL 34",
                                          "MICRO LI 35",
                                          "MICRO LN 35",
                                          "MICRO SP 35",
                                          "RNA LAGU 40",
                                          "RNA BAOU 40",
                                          "RNA NDAM 40",
                                          "RNA BORG 40")),
                x = 1,
                y = 2)

#### set colours
cols <- c(pal[28],
          pal[17],
          pal[9],
          pal[21],
          pal[8],
          pal2[[1]],
          pal[5],
          pal[4])

#### draw legend plot
legend_plot <- ggplot2::ggplot(upset_data,
                               ggplot2::aes(x = x,
                                            y = y,
                                            colour = label,
                                            fill = label)) +
  ggplot2::geom_col(linewidth = 1) +
  ggplot2::scale_colour_manual(name = "Module",
                               values = cols) +
  ggplot2::scale_fill_manual(name = "Module",
                             values = cols) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(colour = ggplot2::alpha(cols,
                                                                                             0.5),
                                                                     stroke = 2))) +
  ggplot2::theme_light()

#### set layout
layout <- "AB
           CE"

#### combine upset plot and legend
plot <- upset_plot +
  legend_plot +
  patchwork::plot_layout(design = layout,
                         guides = "collect") &
  ggplot2::theme(legend.box = "horizontal")

#### save upset plot
png("08_network-analysis/03_figures/04_upset/module-upset.png",
    height = 6,
    res = 300,
    units = "in",
    width = 9)
print(plot)
dev.off()
pdf("08_network-analysis/03_figures/04_upset/module-upset.pdf",
    height = 6,
    width = 9)
print(plot)
dev.off()