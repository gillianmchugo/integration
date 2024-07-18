# combine network enrichment plots
## combine enrichment plots of network results
### install required packages
if (!require(khroma)) install.packages("khroma")
if (!require(magick)) install.packages("magick")
if (!require(magrittr)) install.packages("magrittr")

### function to generate colour palette
smooth_rainbow <- khroma::colour("smooth rainbow")

### generate colour palettes
pal <- smooth_rainbow(31,
                      range = c(0.1,
                                0.9))

pal2 <- khroma::colour("sunset")(11)

### read, crop and add space to legend
legend <- magick::image_read("09_functional-enrichment-analysis/03_network-analysis/03_figures/microarray-modules-01-gprofiler.png") %>%
  magick::image_crop("2700x120")

legend_space <- magick::image_append(c(magick::image_blank(20,
                                                           120,
                                                           color = "white"),
                                       legend))

### read, crop and annotate plots
a <- magick::image_read("09_functional-enrichment-analysis/03_network-analysis/03_figures/microarray-modules-01-gprofiler.png") %>%
  magick::image_crop("2620x1540+80+118")

a_label <- magick::image_append(c(magick::image_blank(100,
                                                      1540,
                                                      color = "white"),
                                  a)) %>%
  magick::image_annotate("Microarray",
                         degrees = 270,
                         gravity = "west",
                         location = "+45+145",
                         size = 60) %>%
  magick::image_composite(magick::image_blank(5,
                                              1450,
                                              color = "#666666"),
                          gravity = "west",
                          offset = "+100+19") %>%
  magick::image_composite(magick::image_blank(5,
                                              486,
                                              color = pal[28]),
                          gravity = "east",
                          offset = "+105-505") %>%
  magick::image_composite(magick::image_blank(5,
                                              486,
                                              color = pal[9]),
                          gravity = "east",
                          offset = "+105+3") %>%
  magick::image_composite(magick::image_blank(5,
                                              486,
                                              color = pal[21]),
                          gravity = "east",
                          offset = "+105+511")

b <- magick::image_read("09_functional-enrichment-analysis/03_network-analysis/03_figures/rna-modules-01-gprofiler.png") %>%
  magick::image_crop("2620x1540+80+118")

b_label <- magick::image_append(c(magick::image_blank(100,
                                                      1540,
                                                      color = "white"),
                                  b)) %>%
  magick::image_annotate("RNA-seq",
                         degrees = 270,
                         gravity = "west",
                         location = "+45+135",
                         size = 60) %>%
  magick::image_composite(magick::image_blank(5,
                                              1450,
                                              color = "#666666"),
                          gravity = "west",
                          offset = "+100+19") %>%
  magick::image_composite(magick::image_blank(5,
                                              359,
                                              color = pal[8]),
                          gravity = "east",
                          offset = "+105-568") %>%
  magick::image_composite(magick::image_blank(5,
                                              359,
                                              color = pal2[[1]]),
                          gravity = "east",
                          offset = "+105-187") %>%
  magick::image_composite(magick::image_blank(5,
                                              359,
                                              color = pal[5]),
                          gravity = "east",
                          offset = "+105+194") %>%
  magick::image_composite(magick::image_blank(5,
                                              359,
                                              color = pal[4]),
                          gravity = "east",
                          offset = "+105+575")

### read, crop and add space to axis titles
x <- magick::image_read("09_functional-enrichment-analysis/03_network-analysis/03_figures/microarray-modules-01-gprofiler.png") %>%
  magick::image_crop("2700x140+0+1660")

x_space <- magick::image_append(c(magick::image_blank(20,
                                                      120,
                                                      color = "white"),
                                  x))

y <- magick::image_read("09_functional-enrichment-analysis/03_network-analysis/03_figures/microarray-modules-01-gprofiler.png") %>%
  magick::image_crop("80x1800")

y_space <- magick::image_composite(magick::image_blank(80,
                                                       3340,
                                                       color = "white"),
                                   y,
                                   gravity = "west")

### combine plots
stack <- magick::image_append(c(legend_space,
                                a_label,
                                b_label,
                                x_space),
                              stack = T)

plot <- magick::image_append(c(y_space,
                               stack))
### save plot
magick::image_write(plot,
                    "09_functional-enrichment-analysis/03_network-analysis/03_figures/gprofiler-modules.png")