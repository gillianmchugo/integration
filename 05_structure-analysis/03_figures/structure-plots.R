# structure plots
## plot structure results
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggh4x)) install.packages("ggh4x")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggtext)) install.packages("ggtext")
if (!require(khroma)) install.packages("khroma")
if (!require(patchwork)) install.packages("patchwork")
if (!require(readr)) install.packages("readr")
if (!require(stringr)) install.packages("stringr")
if (!require(tidyr)) install.packages("tidyr")

### function to read and plot structure data highlighting samples with microarray data
structure_plot_microarray <- function(i, colours, data = "hd", microarray = T, order, tick1, tick2){
  # read and edit fam file
  fam <- readr::read_delim(paste0("01_data-sources/02_filtered-snp-data/",
                                  data,
                                  ".fam"),
                           col_names = F,
                           show_col_types = F) %>%
    dplyr::rename(population = X1,
                  id = X2) %>%
    dplyr::mutate(microarray_data = stringr::str_detect(id,
                                                        "BORA_mchugo|NDAM_mchugo"),
                  group = dplyr::case_when(population %in% c("HOLS", "ANGU", "JERS") ~ "European Bos taurus",
                                           population %in% c("ROMA", "CHIA", "MARC", "MARE", "ALEN") ~ "European hybrid",
                                           population %in% c("MUTU", "LAGU", "NDAG") ~ "African Bos taurus",
                                           population %in% c("NDAM", "BORG", "SOMB", "KETE", "SHEK") ~ "Trypanotolerant African hybrid",
                                           population %in% c("ANKO", "NGAN", "EASZ", "KARA", "BORA") ~ "Trypanosusceptible African hybrid",
                                           population %in% c("THAR", "GIR", "NELO") ~ "Bos indicus")) %>%
    dplyr::mutate(population = factor(population,
                                      levels = c("HOLS", "ANGU", "JERS",
                                                 "ROMA", "CHIA", "MARC", "MARE", "ALEN",
                                                 "MUTU", "LAGU", "NDAG",
                                                 "NDAM", "BORG", "SOMB", "KETE", "SHEK",
                                                 "ANKO", "NGAN", "EASZ", "KARA", "BORA",
                                                 "THAR", "GIR", "NELO")),
                  group = factor(group,
                                 labels = c("European<br>*Bos taurus*",
                                            "European<br>hybrid",
                                            "African<br>*Bos taurus*",
                                            "Trypanotolerant<br>African hybrid",
                                            "Trypanosusceptible<br>African hybrid",
                                            "*Bos<br>indicus*"),
                                 levels = c("European Bos taurus",
                                            "European hybrid",
                                            "African Bos taurus",
                                            "Trypanotolerant African hybrid",
                                            "Trypanosusceptible African hybrid",
                                            "Bos indicus")))
  # set length of custom axis ticks
  ticks <- dplyr::tibble(population = fam$population,
                         id = fam$id,
                         group = fam$group,
                         length = c(rep(NA, 2), tick2, rep(NA, 3),
                                    rep(NA, 12), tick2, rep(NA, 13),
                                    rep(NA, 12), tick1, rep(NA, 12),
                                    rep(NA, 44), tick1, rep(NA, 45),
                                    rep(NA, 24), tick1, rep(NA, 25),
                                    rep(NA, 9), tick1, rep(NA, 9),
                                    rep(NA, 15), tick1, rep(NA, 95),
                                    rep(NA, 13), tick1, rep(NA, 14),
                                    rep(NA, 17), tick1, rep(NA, 18),
                                    rep(NA, 11), tick1, rep(NA, 11),
                                    rep(NA, 7), tick2, rep(NA, 8),
                                    rep(NA, 10), tick1, rep(NA, 11),
                                    rep(NA, 2), tick2, rep(NA, 2),
                                    rep(NA, 6), tick2, rep(NA, 6),
                                    rep(NA, 2), tick1, rep(NA, 2),
                                    rep(NA, 9), tick1, rep(NA, 9),
                                    rep(NA, 23), tick1, rep(NA, 23),
                                    rep(NA, 20), tick2, rep(NA, 20),
                                    rep(NA, 18), tick2, rep(NA, 19),
                                    rep(NA, 13), tick2, rep(NA, 13),
                                    rep(NA, 25), tick2, rep(NA, 25),
                                    rep(NA, 7), tick2, rep(NA, 8),
                                    rep(NA, 11), tick2, rep(NA, 11),
                                    rep(NA, 6), tick2, rep(NA, 6)))
  # read and plot structure results
  readr::read_delim(paste0("01_data-sources/08_structure-results/",
                           dplyr::case_when(data == "hd" ~ "01_hd",
                                            data == "ld" ~ "02_ld"),
                           "/fS_run_K.",
                           i,
                           ".meanQ"),
                    col_names = F,
                    delim = "  ",
                    show_col_types = F) %>%
    dplyr::mutate(population = fam$population,
                  id = fam$id,
                  group = fam$group,
                  microarray_data = fam$microarray_data) %>%
    dplyr::rename(all_of(order)) %>%
    tidyr::pivot_longer(cols = 1:all_of(i)) %>%
    ggplot2::ggplot(ggplot2::aes(x = id,
                                 y = value,
                                 alpha = dplyr::case_when(microarray == T ~ microarray_data,
                                                          microarray == F ~ NA),
                                 fill = name)) +
    ggplot2::geom_bar(stat = "identity",
                      position = "fill",
                      width = 1) +
    ggplot2::scale_alpha_discrete(range = c(0.5,
                                            1)) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggh4x::facet_nested(~ group + population,
                        scales = "free_x",
                        space = "free_x",
                        switch = "x",
                        nest_line = ggplot2::element_line(colour = "grey70",
                                                          linewidth = ggplot2::rel(0.25)),
                        strip = ggh4x::strip_nested(clip = "off",
                                                    text_x = list(ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.6,
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 1,
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.6,
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.8,
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.3,
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = dplyr::case_when(microarray == T ~ "black",
                                                                                                                  microarray == F ~ "grey30"),
                                                                                        face = dplyr::case_when(microarray == T ~ "bold",
                                                                                                                microarray == F ~ "plain"),
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.8,
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.4,
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = dplyr::case_when(microarray == T ~ "black",
                                                                                                                  microarray == F ~ "grey30"),
                                                                                        face = dplyr::case_when(microarray == T ~ "bold",
                                                                                                                microarray == F ~ "plain"),
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.3,
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2)))) +
    ggplot2::geom_linerange(data = ticks,
                            ggplot2::aes(x = id,
                                         ymax = length,
                                         ymin = 0),
                            colour = "grey70",
                            inherit.aes = F,
                            linewidth = rel(0.25),
                            na.rm = T) +
    ggplot2::coord_cartesian(clip = "off",
                             expand = F,
                             ylim = c(0,
                                      1)) +
    ggplot2::ggtitle(paste0("*K* = ",
                            i)) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.border = ggplot2::element_rect(linewidth = ggplot2::rel(0.5)),
                   plot.title = ggtext::element_markdown(size = ggplot2::rel(1)),
                   plot.margin = ggplot2::unit(c(5.5,
                                                 5.5,
                                                 -2.5,
                                                 5.5),
                                               "pt"),
                   strip.background = ggplot2::element_blank())
}

### function to read and plot structure data highlighting samples with expression data
structure_plot_expression <- function(i, colours, expression = T, order, tick1, tick2){
  # read and edit fam file
  fam <- readr::read_delim("03_filter-data/05_filter-snps/ld.fam",
                           col_names = F,
                           show_col_types = F) %>%
    dplyr::rename(population = X1,
                  id = X2) %>%
    dplyr::mutate(expression_data = stringr::str_detect(id,
                                                        "berthier|BORA_mchugo|NDAM_mchugo"),
                  group = dplyr::case_when(population %in% c("HOLS", "ANGU", "JERS") ~ "European Bos taurus",
                                           population %in% c("ROMA", "CHIA", "MARC", "MARE", "ALEN") ~ "European hybrid",
                                           population %in% c("MUTU", "LAGU", "NDAG", "BAOU") ~ "African Bos taurus",
                                           population %in% c("NDAM", "BORG", "SOMB", "KETE", "SHEK") ~ "Trypanotolerant African hybrid",
                                           population %in% c("FULA", "ANKO", "NGAN", "EASZ", "KARA", "BORA") ~ "Trypanosusceptible African hybrid",
                                           population %in% c("THAR", "GIR", "NELO") ~ "Bos indicus")) %>%
    dplyr::mutate(population = factor(population,
                                      levels = c("HOLS", "ANGU", "JERS",
                                                 "ROMA", "CHIA", "MARC", "MARE", "ALEN",
                                                 "MUTU", "LAGU", "NDAG", "BAOU",
                                                 "NDAM", "BORG", "SOMB", "KETE", "SHEK",
                                                 "FULA", "ANKO", "NGAN", "EASZ", "KARA", "BORA",
                                                 "THAR", "GIR", "NELO")),
                  group = factor(group,
                                 labels = c("European<br>*Bos taurus*",
                                            "European<br>hybrid",
                                            "African<br>*Bos taurus*",
                                            "Trypanotolerant<br>African hybrid",
                                            "Trypanosusceptible<br>African hybrid",
                                            "*Bos<br>indicus*"),
                                 levels = c("European Bos taurus",
                                            "European hybrid",
                                            "African Bos taurus",
                                            "Trypanotolerant African hybrid",
                                            "Trypanosusceptible African hybrid",
                                            "Bos indicus")))
  # set length of custom axis ticks
  ticks <- dplyr::tibble(population = fam$population,
                         id = fam$id,
                         group = fam$group,
                         length = c(rep(NA, 2), tick2, rep(NA, 3),
                                    rep(NA, 12), tick2, rep(NA, 13),
                                    rep(NA, 12), tick1, rep(NA, 12),
                                    rep(NA, 18), tick2, rep(NA, 18),
                                    rep(NA, 55), tick1, rep(NA, 54),
                                    rep(NA, 130), tick2, rep(NA, 130),
                                    rep(NA, 9), tick1, rep(NA, 9),
                                    rep(NA, 15), tick1, rep(NA, 95),
                                    rep(NA, 25), tick2, rep(NA, 25),
                                    rep(NA, 13), tick1, rep(NA, 14),
                                    rep(NA, 17), tick1, rep(NA, 18),
                                    rep(NA, 11), tick1, rep(NA, 11),
                                    rep(NA, 7), tick2, rep(NA, 8),
                                    rep(NA, 10), tick2, rep(NA, 11),
                                    rep(NA, 25), tick2, rep(NA, 25),
                                    rep(NA, 6), tick2, rep(NA, 6),
                                    rep(NA, 2), tick1, rep(NA, 2),
                                    rep(NA, 9), tick1, rep(NA, 9),
                                    rep(NA, 23), tick1, rep(NA, 23),
                                    rep(NA, 40), tick1, rep(NA, 39),
                                    rep(NA, 18), tick2, rep(NA, 19),
                                    rep(NA, 13), tick2, rep(NA, 13),
                                    rep(NA, 25), tick2, rep(NA, 25),
                                    rep(NA, 7), tick1, rep(NA, 8),
                                    rep(NA, 11), tick1, rep(NA, 11),
                                    rep(NA, 6), tick2, rep(NA, 6)))
  # read and plot structure results
  readr::read_delim(paste0("05_structure-analysis/02_output/fS_run_K.",
                           i,
                           ".meanQ"),
                    col_names = F,
                    delim = "  ",
                    show_col_types = F) %>%
    dplyr::mutate(population = fam$population,
                  id = fam$id,
                  group = fam$group,
                  expression_data = fam$expression_data) %>%
    dplyr::rename(all_of(order)) %>%
    tidyr::pivot_longer(cols = 1:all_of(i)) %>%
    ggplot2::ggplot(ggplot2::aes(x = id,
                                 y = value,
                                 alpha = dplyr::case_when(expression == T ~ expression_data,
                                                          expression == F ~ NA),
                                 fill = name)) +
    ggplot2::geom_bar(stat = "identity",
                      position = "fill",
                      width = 1) +
    ggplot2::scale_alpha_discrete(range = c(0.5,
                                            1)) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_y_continuous(expand = c(0,
                                           0)) +
    ggh4x::facet_nested(~ group + population,
                        scales = "free_x",
                        space = "free_x",
                        switch = "x",
                        nest_line = ggplot2::element_line(colour = "grey70",
                                                          linewidth = ggplot2::rel(0.25)),
                        strip = ggh4x::strip_nested(clip = "off",
                                                    text_x = list(ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 1,
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.8,
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.8,
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.6,
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.4,
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.1,
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = dplyr::case_when(expression == T ~ "black",
                                                                                                                  expression == F ~ "grey30"),
                                                                                        face = dplyr::case_when(expression == T ~ "bold",
                                                                                                                expression == F ~ "plain"),
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = dplyr::case_when(expression == T ~ "black",
                                                                                                                  expression == F ~ "grey30"),
                                                                                        face = dplyr::case_when(expression == T ~ "bold",
                                                                                                                expression == F ~ "plain"),
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = dplyr::case_when(expression == T ~ "black",
                                                                                                                  expression == F ~ "grey30"),
                                                                                        face = dplyr::case_when(expression == T ~ "bold",
                                                                                                                expression == F ~ "plain"),
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = dplyr::case_when(expression == T ~ "black",
                                                                                                                  expression == F ~ "grey30"),
                                                                                        face = dplyr::case_when(expression == T ~ "bold",
                                                                                                                expression == F ~ "plain"),
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = dplyr::case_when(expression == T ~ "black",
                                                                                                                  expression == F ~ "grey30"),
                                                                                        face = dplyr::case_when(expression == T ~ "bold",
                                                                                                                expression == F ~ "plain"),
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = dplyr::case_when(expression == T ~ "black",
                                                                                                                  expression == F ~ "grey30"),
                                                                                        face = dplyr::case_when(expression == T ~ "bold",
                                                                                                                expression == F ~ "plain"),
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2)))) +
    ggplot2::geom_linerange(data = ticks,
                            ggplot2::aes(x = id,
                                         ymax = length,
                                         ymin = 0),
                            colour = "grey70",
                            inherit.aes = F,
                            linewidth = ggplot2::rel(0.25),
                            na.rm = T) +
    ggplot2::coord_cartesian(clip = "off",
                             expand = F,
                             ylim = c(0,
                                      1)) +
    ggplot2::ggtitle(paste0("*K* = ",
                            i)) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.border = ggplot2::element_rect(linewidth = ggplot2::rel(0.5)),
                   plot.title = ggtext::element_markdown(size = ggplot2::rel(1)),
                   plot.margin = ggplot2::unit(c(5.5,
                                                 5.5,
                                                 -2.5,
                                                 5.5),
                                               "pt"),
                   strip.background = ggplot2::element_blank())
}

### function to save structure plot for individual k value highlighting samples with microarray data
individual_k_plot_microarray <- function(i, colours, data = "hd", order, tick1 = -0.01, tick2 = -0.04){
  plot <- structure_plot_microarray(i,
                                    colours = colours,
                                    data = data,
                                    order = order,
                                    tick1 = tick1,
                                    tick2 = tick2)
  png(paste0("05_structure-analysis/03_figures/",
             data,
             "-k-",
             stringr::str_pad(i,
                              2,
                              pad = "0",
                              side = "left"),
             "-microarray.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("05_structure-analysis/03_figures/",
             data,
             "-k-",
             stringr::str_pad(i,
                              2,
                              pad = "0",
                              side = "left"),
             "-microarray.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
  plot <- structure_plot_microarray(i,
                                    colours = colours,
                                    data = data,
                                    microarray = F,
                                    order = order,
                                    tick1 = tick1,
                                    tick2 = tick2)
  png(paste0("05_structure-analysis/03_figures/",
             data,
             "-k-",
             stringr::str_pad(i,
                              2,
                              pad = "0",
                              side = "left"),
             ".png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("05_structure-analysis/03_figures/",
             data,
             "-k-",
             stringr::str_pad(i,
                              2,
                              pad = "0",
                              side = "left"),
             ".pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

### function to save structure plot for individual k value highlighting samples with expression data
individual_k_plot_expression <- function(i, colours, order, tick1 = -0.01, tick2 = -0.04){
  plot <- structure_plot_expression(i,
                                    colours = colours,
                                    order = order,
                                    tick1 = tick1,
                                    tick2 = tick2)
  png(paste0("05_structure-analysis/03_figures/",
             "k-",
             stringr::str_pad(i,
                              2,
                              pad = "0",
                              side = "left"),
             "-expression.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("05_structure-analysis/03_figures/",
             "k-",
             stringr::str_pad(i,
                              2,
                              pad = "0",
                              side = "left"),
             "-expression.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
  plot <- structure_plot_expression(i,
                                    colours = colours,
                                    expression = F,
                                    order = order,
                                    tick1 = tick1,
                                    tick2 = tick2)
  png(paste0("05_structure-analysis/03_figures/",
             "k-",
             stringr::str_pad(i,
                              2,
                              pad = "0",
                              side = "left"),
             ".png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("05_structure-analysis/03_figures/",
             "k-",
             stringr::str_pad(i,
                              2,
                              pad = "0",
                              side = "left"),
             ".pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

### function to generate colour palette
smooth_rainbow <- khroma::colour("smooth rainbow")

### generate colour palette
pal <- smooth_rainbow(31,
                      range = c(0.1,
                                0.9))

### assign colours and orders
hd_02_colours <- c(pal[17], pal[28])
hd_02_order <- c(a = "X2", b = "X1")

ld_02_colours <- c(pal[17], pal[28])
ld_02_order <- c(a = "X2", b = "X1")

k_02_colours <- c(pal[9], pal[28])
k_02_order <- c(a = "X2", b = "X1")

hd_03_colours <- c(pal[9], pal[17], pal[28])
hd_03_order <- c(a = "X1", b = "X3", c = "X2")

ld_03_colours <- c(pal[9], pal[17], pal[28])
ld_03_order <- c(a = "X1", b = "X2", c = "X3")

k_03_colours <- c(pal[9], pal[17], pal[28])
k_03_order <- c(a = "X1", b = "X2", c = "X3")

### apply functions to results
#### k = 2
individual_k_plot_microarray(2,
                             colours = hd_02_colours,
                             order = hd_02_order)

individual_k_plot_microarray(2,
                             colours = ld_02_colours,
                             data = "ld",
                             order = ld_02_order)

individual_k_plot_expression(2,
                             colours = k_02_colours,
                             order = k_02_order)

#### k = 3
individual_k_plot_microarray(3,
                             colours = hd_03_colours,
                             order = hd_03_order)

individual_k_plot_microarray(3,
                             colours = ld_03_colours,
                             data = "ld",
                             order = ld_03_order)

individual_k_plot_expression(3,
                             colours = k_03_colours,
                             order = k_03_order)