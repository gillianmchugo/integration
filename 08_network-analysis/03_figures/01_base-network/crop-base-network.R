# crop base network
## crop base network plot
### install required packages
if (!require(magick)) install.packages("magick")
if (!require(magrittr)) install.packages("magrittr")

### read and crop base network
base_network <- magick::image_read_pdf("08_network-analysis/03_figures/01_base-network/base-network.pdf") %>%
  magick::image_crop("2300x2200+604")

### save plot
magick::image_write(base_network,
                    "08_network-analysis/03_figures/01_base-network/base-network-crop.png")