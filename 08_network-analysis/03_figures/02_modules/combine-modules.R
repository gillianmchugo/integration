# combine modules
## combine functional modules from network analysis
### install required packages
if (!require(magick)) install.packages("magick")
if (!require(magrittr)) install.packages("magrittr")
if (!require(pdftools)) install.packages("pdftools")

### read legend
legend <- magick::image_read("08_network-analysis/03_figures/02_modules/legend.png") %>%
  magick::image_scale("9200")

### read and label modules a-d
a <- magick::image_read_pdf("08_network-analysis/03_figures/02_modules/MICRO-BL-34-module-01.pdf") %>%
  magick::image_crop("2300x2200+604")

a_label <- magick::image_append(c(magick::image_blank(2300,
                                                      140,
                                                      color = "white"),
                                  a),
                                stack = T) %>%
  magick::image_annotate("MICRO BL 34",
                         gravity = "north",
                         size = 140)

b <- magick::image_read_pdf("08_network-analysis/03_figures/02_modules/MICRO-LI-35-module-01.pdf") %>%
  magick::image_crop("2300x2200+604")

b_label <- magick::image_append(c(magick::image_blank(2300,
                                                      140,
                                                      color = "white"),
                                  b),
                                stack = T) %>%
  magick::image_annotate("MICRO LI 35",
                         gravity = "north",
                         size = 140)

c <- magick::image_read_pdf("08_network-analysis/03_figures/02_modules/RNA-LAGU-40-module-01.pdf") %>%
  magick::image_crop("2300x2200+604")

c_label <- magick::image_append(c(magick::image_blank(2300,
                                                      140,
                                                      color = "white"),
                                  c),
                                stack = T) %>%
  magick::image_annotate("RNA LAGU 40",
                         gravity = "north",
                         size = 140)

d <- magick::image_read_pdf("08_network-analysis/03_figures/02_modules/RNA-BAOU-40-module-01.pdf") %>%
  magick::image_crop("2300x2200+604")

d_label <- magick::image_append(c(magick::image_blank(2300,
                                                      140,
                                                      color = "white"),
                                  d),
                                stack = T) %>%
  magick::image_annotate("RNA BAOU 40",
                         gravity = "north",
                         size = 140)

### combine modules a-d
abcd <- magick::image_append(c(a_label,
                               b_label,
                               c_label,
                               d_label))

### read and label modules e-h
e <- magick::image_read_pdf("08_network-analysis/03_figures/02_modules/MICRO-LN-35-module-01.pdf") %>%
  magick::image_crop("2300x2200+604")

e_label <- magick::image_append(c(magick::image_blank(2300,
                                                      140,
                                                      color = "white"),
                                  e),
                                stack = T) %>%
  magick::image_annotate("MICRO LN 35",
                         gravity = "north",
                         size = 140)

f <- magick::image_read_pdf("08_network-analysis/03_figures/02_modules/MICRO-SP-35-module-01.pdf") %>%
  magick::image_crop("2300x2200+604")

f_label <- magick::image_append(c(magick::image_blank(2300,
                                                      140,
                                                      color = "white"),
                                  f),
                                stack = T) %>%
  magick::image_annotate("MICRO SP 35",
                         gravity = "north",
                         size = 140)

ef <- magick::image_append(c(e_label,
                             f_label))

ef_label <- magick::image_append(c(ef,
                                   magick::image_blank(4600,
                                                       190,
                                                       color = "white")),
                                 stack = T) %>%
  magick::image_annotate("Micorarray",
                         gravity = "south",
                         size = 140) %>%
  magick::image_composite(magick::image_blank(4000,
                                              10,
                                              color = "#666666"),
                          gravity = "south",
                          offset = "+0+180")

g <- magick::image_read_pdf("08_network-analysis/03_figures/02_modules/RNA-NDAM-40-module-01.pdf") %>%
  magick::image_crop("2300x2200+604")

g_label <- magick::image_append(c(magick::image_blank(2300,
                                                      140,
                                                      color = "white"),
                                  g),
                                stack = T) %>%
  magick::image_annotate("RNA NDAM 40",
                         gravity = "north",
                         size = 140)

h <- magick::image_read_pdf("08_network-analysis/03_figures/02_modules/RNA-BORG-40-module-01.pdf") %>%
  magick::image_crop("2300x2200+604")

h_label <- magick::image_append(c(magick::image_blank(2300,
                                                      140,
                                                      color = "white"),
                                  h),
                                stack = T) %>%
  magick::image_annotate("RNA BORG 40",
                         gravity = "north",
                         size = 140)

gh <- magick::image_append(c(g_label,
                             h_label))

gh_label <- magick::image_append(c(gh,
                                   magick::image_blank(4600,
                                                       190,
                                                       color = "white")),
                                 stack = T) %>%
  magick::image_annotate("RNA-seq",
                         gravity = "south",
                         size = 140) %>%
  magick::image_composite(magick::image_blank(4000,
                                              10,
                                              color = "#666666"),
                          gravity = "south",
                          offset = "+0+180")

### combine modules e-h
efgh <- magick::image_append(c(ef_label,
                               gh_label))

### combine all modules and legend
plot <- magick::image_append(c(legend,
                               abcd,
                               efgh),
                             stack = T)

### save plot
magick::image_write(plot,
                    "08_network-analysis/03_figures/02_modules/all-modules.png")