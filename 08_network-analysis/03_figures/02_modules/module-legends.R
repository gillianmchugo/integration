# module legends
## add legend to module plots
### install required packages
if (!require(magick)) install.packages("magick")
if (!require(magrittr)) install.packages("magrittr")

### function to add legend to module plots
module_legend <- function(module){
  # read legend
  legend <- magick::image_read("08_network-analysis/03_figures/02_modules/legend.png") %>%
    magick::image_scale("2300")
  # read module
  network <- magick::image_read_pdf(paste0("08_network-analysis/03_figures/02_modules/",
                                           module,
                                           "-module-01.pdf")) %>%
    magick::image_crop("2300x2200+604")
  # combine module and legend
  plot <- magick::image_append(c(legend,
                                 network),
                               stack = T)
  # save plot
  magick::image_write(plot, 
                      paste0("08_network-analysis/03_figures/02_modules/",
                             module,
                             "-module-01-legend.png"))
}

### set modules
modules <- c("MICRO-BL-34",
             "MICRO-LI-35",
             "MICRO-LN-35",
             "MICRO-SP-35",
             "RNA-LAGU-40",
             "RNA-BAOU-40",
             "RNA-NDAM-40",
             "RNA-BORG-40")

### apply function
lapply(modules,
       module_legend)