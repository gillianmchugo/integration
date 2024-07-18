# combine local ancestry plots
## combine local ancestry plots to compare analyses
### install required packages
if (!require(khroma)) install.packages("khroma")
if (!require(magick)) install.packages("magick")
if (!require(magrittr)) install.packages("magrittr")

### function to combine local ancestry plots
combine_local_ancestry_plots <- function(expression = F){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palette
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # read legend
  legend <- magick::image_read(paste0("06_local-ancestry-analysis/03_figures/01_BAOU/elai-ld-BAOU-01",
                                      dplyr::case_when(expression == T ~ "-expression",
                                                       expression == F ~ ""),
                                      ".png")) %>%
    magick::image_crop("2700x170")
  # read and label local ancestry plots a-c
  a <- magick::image_read(paste0("06_local-ancestry-analysis/03_figures/05_LAGU/elai-ld-LAGU",
                                 dplyr::case_when(expression == T ~ "-expression",
                                                  expression == F ~ ""),
                                 "-round.png"))%>%
    magick::image_crop("1800x1800+200") %>%
    magick::image_annotate("LAGU",
                           gravity = "center",
                           size = 100)
  b <- magick::image_read(paste0("06_local-ancestry-analysis/03_figures/06_NDAM/elai-ld-NDAM",
                                 dplyr::case_when(expression == T ~ "-expression",
                                                  expression == F ~ ""),
                                 "-round.png")) %>%
    magick::image_crop("1800x1800+200") %>%
    magick::image_annotate("NDAM",
                           gravity = "center",
                           size = 100)
  c <- magick::image_read(paste0("06_local-ancestry-analysis/03_figures/04_FULA/elai-ld-FULA",
                                 dplyr::case_when(expression == T ~ "-expression",
                                                  expression == F ~ ""),
                                 "-round.png")) %>%
    magick::image_crop("1800x1800+200") %>%
    magick::image_annotate("FULA",
                           gravity = "center",
                           size = 100)
  # combine local ancestry plots a-c
  abc <- magick::image_append(c(a,
                                b,
                                c)) %>%
    magick::image_scale("2700")
  # read and label local ancestry plots d-f
  d <- magick::image_read(paste0("06_local-ancestry-analysis/03_figures/01_BAOU/elai-ld-BAOU",
                                 dplyr::case_when(expression == T ~ "-expression",
                                                  expression == F ~ ""),
                                 "-round.png")) %>%
    magick::image_crop("1800x1800+200") %>%
    magick::image_annotate("BAOU",
                           gravity = "center",
                           size = 100) %>%
    magick::image_annotate("African",
                           gravity = "south",
                           location = "-255+0",
                           size = 100) %>%
    magick::image_annotate("Bos taurus",
                           gravity = "south",
                           location = "+170+0",
                           size = 100,
                           style = "italic") %>%
    magick::image_composite(magick::image_blank(1400,
                                                10,
                                                color = pal[9]),
                            gravity = "south",
                            offset = "+0+125")
  e <- magick::image_read(paste0("06_local-ancestry-analysis/03_figures/03_BORG/elai-ld-BORG",
                                 dplyr::case_when(expression == T ~ "-expression",
                                                  expression == F ~ ""),
                                 "-round.png")) %>%
    magick::image_crop("1800x1800+200") %>%
    magick::image_annotate("BORG",
                           gravity = "center",
                           size = 100) %>%
    magick::image_annotate("Trypanotolerant African hybrid",
                           gravity = "south",
                           size = 100) %>%
    magick::image_composite(magick::image_blank(1400,
                                                10,
                                                color = pal[5]),
                            gravity = "south",
                            offset = "+0+125")
  f <- magick::image_read(paste0("06_local-ancestry-analysis/03_figures/02_BORA/elai-ld-BORA",
                                 dplyr::case_when(expression == T ~ "-expression",
                                                  expression == F ~ ""),
                                 "-round.png")) %>%
    magick::image_crop("1800x1800+200") %>%
    magick::image_annotate("BORA",
                           gravity = "center",
                           size = 100) %>%
    magick::image_annotate("Trypanosusceptible African hybrid",
                           gravity = "south",
                           size = 100) %>%
    magick::image_composite(magick::image_blank(1400,
                                                10,
                                                color = pal[26]),
                            gravity = "south",
                            offset = "+0+125")
  # combine local ancestry plots d-f
  def <- magick::image_append(c(d,
                                e,
                                f)) %>%
    magick::image_scale("2700")
  # combine all plots and legend
  plot <- magick::image_append(c(legend,
                                 abc,
                                 def),
                               stack = T)
  # save plot
  magick::image_write(plot,
                      paste0("06_local-ancestry-analysis/03_figures/07_comparison/local-ancestry-plots-comparison",
                             dplyr::case_when(expression == T ~ "-expression",
                                              expression == F ~ ""),
                             ".png"))
}

### apply function
combine_local_ancestry_plots_expression()
combine_local_ancestry_plots_expression(expression = T)