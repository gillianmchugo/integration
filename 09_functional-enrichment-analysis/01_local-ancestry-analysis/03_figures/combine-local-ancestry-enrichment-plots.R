# combine local ancestry enrichment plots
## combine enrichment plots of local ancestry results
### install required packages
if (!require(khroma)) install.packages("khroma")
if (!require(magick)) install.packages("magick")
if (!require(magrittr)) install.packages("magrittr")

### function to combine gprofiler enrichment plots of local ancestry results
combine_local_ancestry_gprofiler_plots <- function(group, colours, expression = F){
  populations <- get(group)
  cols <- get(paste0(group,
                     "_colours"))
  # read, crop and annotate plots
  a <- magick::image_read(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/03_figures/gprofiler-elai-ld-",
                                 populations[1],
                                 dplyr::case_when(expression == T ~ "-expression",
                                                  expression == F ~ ""),
                                 ".png")) %>%
    magick::image_crop("2620x1540+80+118")
  a_label <- magick::image_append(c(magick::image_blank(100,
                                                        1540,
                                                        color = "white"),
                                    a)) %>%
    magick::image_annotate(populations[1],
                           degrees = 270,
                           gravity = "west",
                           location = "+45+75",
                           size = 60) %>%
    magick::image_composite(magick::image_blank(5,
                                                1450,
                                                color = cols[1]),
                            gravity = "west",
                            offset = "+100+19") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[17]),
                            gravity = "east",
                            offset = "+110-505") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[9]),
                            gravity = "east",
                            offset = "+110+3") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[28]),
                            gravity = "east",
                            offset = "+110+511")
  b <- magick::image_read(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/03_figures/gprofiler-elai-ld-",
                                 populations[2],
                                 dplyr::case_when(expression == T ~ "-expression",
                                                  expression == F ~ ""),
                                 ".png")) %>%
    magick::image_crop("2620x1540+80+118")
  b_label <- magick::image_append(c(magick::image_blank(100,
                                                        1540,
                                                        color = "white"),
                                    b)) %>%
    magick::image_annotate(populations[2],
                           degrees = 270,
                           gravity = "west",
                           location = "+45+75",
                           size = 60) %>%
    magick::image_composite(magick::image_blank(5,
                                                1450,
                                                color = cols[2]),
                            gravity = "west",
                            offset = "+100+19") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[17]),
                            gravity = "east",
                            offset = "+110-505") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[9]),
                            gravity = "east",
                            offset = "+110+3") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[28]),
                            gravity = "east",
                            offset = "+110+511")
  # read, crop and add space to legend and axis titles
  legend <- magick::image_read(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/03_figures/gprofiler-elai-ld-",
                                      populations[1],
                                      dplyr::case_when(expression == T ~ "-expression",
                                                       expression == F ~ ""),
                                      ".png")) %>%
    magick::image_crop("2700x120")
  legend_space <- magick::image_append(c(magick::image_blank(20,
                                                             120,
                                                             color = "white"),
                                         legend))
  x <- magick::image_read(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/03_figures/gprofiler-elai-ld-",
                                 populations[1],
                                 dplyr::case_when(expression == T ~ "-expression",
                                                  expression == F ~ ""),
                                 ".png")) %>%
    magick::image_crop("2700x140+0+1660")
  x_space <- magick::image_append(c(magick::image_blank(20,
                                                        120,
                                                        color = "white"),
                                    x))
  y <- magick::image_read(paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/03_figures/gprofiler-elai-ld-",
                                 populations[1],
                                 dplyr::case_when(expression == T ~ "-expression",
                                                  expression == F ~ ""),
                                 ".png")) %>%
    magick::image_crop("80x1800")
  y_space <- magick::image_composite(magick::image_blank(80,
                                                         3340,
                                                         color = "white"),
                                     y,
                                     gravity = "west")
  # combine plots
  stack <- magick::image_append(c(legend_space,
                                  a_label,
                                  b_label,
                                  x_space),
                                stack = T)
  plot <- magick::image_append(c(y_space,
                                 stack))
  # save plot
  magick::image_write(plot,
                      paste0("09_functional-enrichment-analysis/01_local-ancestry-analysis/03_figures/gprofiler-elai-ld-",
                             stringr::str_replace_all(group,
                                                      "_",
                                                      "-"),
                             dplyr::case_when(expression == T ~ "-expression",
                                              expression == F ~ ""),
                             ".png"))
}

### function to generate colour palette
smooth_rainbow <- khroma::colour("smooth rainbow")

### generate colour palettes
pal <- smooth_rainbow(31,
                      range = c(0.1,
                                0.9))

pal2 <- khroma::colour("sunset")(11)

### set populations and colours
african_bos_taurus <- c("LAGU",
                        "BAOU")

african_bos_taurus_colours <- c(pal[8],
                                pal2[[1]])

trypanotolerant_african_hybrids <- c("NDAM",
                                     "BORG")

trypanotolerant_african_hybrids_colours <- c(pal[5],
                                             pal[4])

trypanosusceptible_african_hybrids <- c("FULA",
                                        "BORA")

trypanosusceptible_african_hybrids_colours <- c(pal[21],
                                                pal[26])

### apply function
combine_local_ancestry_gprofiler_plots("african_bos_taurus")
combine_local_ancestry_gprofiler_plots("african_bos_taurus",
                                       expression = T)

combine_local_ancestry_gprofiler_plots("trypanotolerant_african_hybrids")
combine_local_ancestry_gprofiler_plots("trypanotolerant_african_hybrids",
                                       expression = T)

combine_local_ancestry_gprofiler_plots("trypanosusceptible_african_hybrids")
combine_local_ancestry_gprofiler_plots("trypanosusceptible_african_hybrids",
                                       expression = T)