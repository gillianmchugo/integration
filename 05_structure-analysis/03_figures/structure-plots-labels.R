# structure plots labels
## add labels to plots of structure results
### install required packages
if (!require(khroma)) install.packages("khroma")
if (!require(magick)) install.packages("magick")
if (!require(magrittr)) install.packages("magrittr")

### function to add labels to plots of samples with microarray data
microarray_labels <- function(input){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palette
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # read image and add label annotation
  magick::image_read(paste0("05_structure-analysis/03_figures/",
                            input,
                            ".png")) %>%
    magick::image_composite(magick::image_blank(280,
                                                5,
                                                color = pal[17]),
                            gravity = "southwest",
                            offset = "+110+120") %>%
    magick::image_composite(magick::image_blank(350,
                                                5,
                                                color = pal[15]),
                            gravity = "southwest",
                            offset = "+410+120") %>%
    magick::image_composite(magick::image_blank(245,
                                                5,
                                                color = pal[9]),
                            gravity = "southwest",
                            offset = "+780+120") %>%
    magick::image_composite(magick::image_blank(510,
                                                5,
                                                color = pal[5]),
                            gravity = "south",
                            offset = "-55+120") %>%
    magick::image_composite(magick::image_blank(825,
                                                5,
                                                color = pal[26]),
                            gravity = "southeast",
                            offset = "+305+120") %>%
    magick::image_composite(magick::image_blank(265,
                                                5,
                                                color = pal[28]),
                            gravity = "southeast",
                            offset = "+20+120") %>%
    magick::image_write(paste0("05_structure-analysis/03_figures/",
                               input,
                               "-labels.png"))
}

### function to add labels to plots of samples with expression data
expression_labels <- function(input){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palette
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # read image and add label annotation
  magick::image_read(paste0("05_structure-analysis/03_figures/",
                            input,
                            ".png")) %>%
    magick::image_composite(magick::image_blank(200,
                                                5,
                                                color = pal[17]),
                            gravity = "southwest",
                            offset = "+110+120") %>%
    magick::image_composite(magick::image_blank(260,
                                                5,
                                                color = pal[15]),
                            gravity = "southwest",
                            offset = "+325+120") %>%
    magick::image_composite(magick::image_blank(335,
                                                5,
                                                color = pal[9]),
                            gravity = "southwest",
                            offset = "+605+120") %>%
    magick::image_composite(magick::image_blank(790,
                                                5,
                                                color = pal[5]),
                            gravity = "south",
                            offset = "+5+120") %>%
    magick::image_composite(magick::image_blank(705,
                                                5,
                                                color = pal[26]),
                            gravity = "southeast",
                            offset = "+225+120") %>%
    magick::image_composite(magick::image_blank(185,
                                                5,
                                                color = pal[28]),
                            gravity = "southeast",
                            offset = "+20+120") %>%
    magick::image_write(paste0("05_structure-analysis/03_figures/",
                               input,
                               "-labels.png"))
}

### set microarray inputs
microarray_inputs <- c("hd-k-02",
                       "hd-k-02-microarray",
                       "hd-k-03",
                       "hd-k-03-microarray",
                       "ld-k-02",
                       "ld-k-02-microarray",
                       "ld-k-03",
                       "ld-k-03-microarray")

### set expression inputs
expression_inputs <- c("k-02",
                       "k-02-expression",
                       "k-03",
                       "k-03-expression")

### apply functions
lapply(microarray_inputs,
       microarray_labels)

lapply(expression_inputs,
       expression_labels)