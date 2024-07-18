# local ancestry plots
## plot local ancestry results
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(khroma)) install.packages("khroma")
if (!require(parallel)) install.packages("parallel")
if (!require(patchwork)) install.packages("patchwork")
if (!require(readr)) install.packages("readr")
if (!require(scales)) install.packages("scales")
if (!require(stringr)) install.packages("stringr")
if (!require(tibble)) install.packages("tibble")

### function to plot individual chromosomes for samples with and without expression data
chromosome_plot_expression <- function(chromosome, population, expression = T){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palette
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # plot results
  plot <- readr::read_csv(paste0("06_local-ancestry-analysis/",
                                 "02_output/",
                                 dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                                                  population == "BORA" ~ "02_BORA",
                                                  population == "BORG" ~ "03_BORG",
                                                  population == "FULA" ~ "04_FULA",
                                                  population == "LAGU" ~ "05_LAGU",
                                                  population == "NDAM" ~ "06_NDAM"),
                                 "/elai-ld-",
                                 population,
                                 "-",
                                 stringr::str_pad(chromosome,
                                                  2,
                                                  pad = "0",
                                                  side = "left"),
                                 dplyr::case_when(expression == T ~ "-expression",
                                                  expression == F ~ ""),
                                 ".csv")) %>%
    ggplot2::ggplot(ggplot2::aes(x = pos)) +
    ggplot2::geom_area(ggplot2::aes(y = mean_a +
                                      mean_b +
                                      mean_c,
                                    fill = "a")) +
    ggplot2::geom_area(ggplot2::aes(y = mean_b +
                                      mean_c,
                                    fill = "b")) +
    ggplot2::geom_area(ggplot2::aes(y = mean_c,
                                    fill = "c")) +
    ggplot2::scale_fill_manual(labels = c(expression(paste("European ",
                                                           italic("Bos taurus"),
                                                           phantom("p"))),
                                          expression(paste("African ",
                                                           italic("Bos taurus"),
                                                           phantom("p"))),
                                          expression(paste(italic("Bos indicus"),
                                                           phantom("p")))),
                               values = c(pal[17],
                                          pal[9],
                                          pal[28])) +
    ggplot2::scale_x_continuous(expand=c(0,
                                         0),
                                labels = scales::label_number(scale_cut = c(0,
                                                                            "Mb" = 1000000)),
                                name = paste("Position on Chromosome",
                                             chromosome)) +
    ggplot2::scale_y_continuous(expand=c(0,
                                         0),
                                name = "Mean Ancestry") +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "top",
                   legend.text.align = 0,
                   legend.title = ggplot2::element_blank())
  # save plot
  png(paste0("06_local-ancestry-analysis/",
             "03_figures/",
             dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                              population == "BORA" ~ "02_BORA",
                              population == "BORG" ~ "03_BORG",
                              population == "FULA" ~ "04_FULA",
                              population == "LAGU" ~ "05_LAGU",
                              population == "NDAM" ~ "06_NDAM"),
             "/elai-ld-",
             population,
             "-",
             stringr::str_pad(chromosome,
                              2,
                              pad = "0",
                              side = "left"),
             dplyr::case_when(expression == T ~ "-expression",
                              expression == F ~ ""),
             ".png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("06_local-ancestry-analysis/",
             "03_figures/",
             dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                              population == "BORA" ~ "02_BORA",
                              population == "BORG" ~ "03_BORG",
                              population == "FULA" ~ "04_FULA",
                              population == "LAGU" ~ "05_LAGU",
                              population == "NDAM" ~ "06_NDAM"),
             "/elai-ld-",
             population,
             "-",
             stringr::str_pad(chromosome,
                              2,
                              pad = "0",
                              side = "left"),
             dplyr::case_when(expression == T ~ "-expression",
                              expression == F ~ ""),
             ".pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

### function to plot flat genome data for samples with and without expression data
flat_genome_plot_expression <- function(population, expression = T){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palette
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # get input data
  input <- readr::read_csv(paste0("06_local-ancestry-analysis/",
                                  "02_output/",
                                  dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                                                   population == "BORA" ~ "02_BORA",
                                                   population == "BORG" ~ "03_BORG",
                                                   population == "FULA" ~ "04_FULA",
                                                   population == "LAGU" ~ "05_LAGU",
                                                   population == "NDAM" ~ "06_NDAM"),
                                  "/elai-ld-",
                                  population,
                                  dplyr::case_when(expression == T ~ "-expression",
                                                   expression == F ~ ""),
                                  ".csv"))
  # get chromosome lengths
  length <- input %>%
    dplyr::count(chr)
  # set label positions
  label_pos <- c()
  for(i in 1:29){
    label_pos <- c(label_pos,
                   sum(length$n[0:(i-1)],
                       length$n[i]/2))
  }
  # plot results
  plot <- ggplot2::ggplot(input,
                          ggplot2::aes(x = row_num)) +
    ggplot2::geom_area(ggplot2::aes(y = mean_a +
                                      mean_b +
                                      mean_c,
                                    fill = "a")) +
    ggplot2::geom_area(ggplot2::aes(y = mean_b +
                                      mean_c,
                                    fill = "b")) +
    ggplot2::geom_area(ggplot2::aes(y = mean_c,
                                    fill = "c")) +
    ggplot2::scale_fill_manual(labels = c(expression(paste("European ",
                                                           italic("Bos taurus"),
                                                           phantom("p"))),
                                          expression(paste("African ",
                                                           italic("Bos taurus"),
                                                           phantom("p"))),
                                          expression(paste(italic("Bos indicus"),
                                                           phantom("p")))),
                               values = c(pal[17],
                                          pal[9],
                                          pal[28])) +
    ggplot2::scale_x_continuous(breaks = label_pos,
                                expand = c(0,
                                           0),
                                labels = 1:29,
                                name = "Chromosome") +
    ggplot2::scale_y_continuous(expand = c(0,
                                           0),
                                name = "Mean Ancestry") +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "top",
                   legend.text.align = 0,
                   legend.title = ggplot2::element_blank())
  # add line segments between chromosomes
  for(i in 1:28){
    plot <- plot +
      ggplot2::geom_segment(x = sum(length$n[1:i],
                                    1),
                            xend = sum(length$n[1:i],
                                       1),
                            y = 0,
                            yend = 1,
                            colour = "grey70",
                            linewidth = ggplot2::rel(0.25))
  }
  # save plot
  png(paste0("06_local-ancestry-analysis/",
             "03_figures/",
             dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                              population == "BORA" ~ "02_BORA",
                              population == "BORG" ~ "03_BORG",
                              population == "FULA" ~ "04_FULA",
                              population == "LAGU" ~ "05_LAGU",
                              population == "NDAM" ~ "06_NDAM"),
             "/elai-ld-",
             population,
             dplyr::case_when(expression == T ~ "-expression",
                              expression == F ~ ""),
             "-flat.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("06_local-ancestry-analysis/",
             "03_figures/",
             dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                              population == "BORA" ~ "02_BORA",
                              population == "BORG" ~ "03_BORG",
                              population == "FULA" ~ "04_FULA",
                              population == "LAGU" ~ "05_LAGU",
                              population == "NDAM" ~ "06_NDAM"),
             "/elai-ld-",
             population,
             dplyr::case_when(expression == T ~ "-expression",
                              expression == F ~ ""),
             "-flat.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

### function to plot round genome data for samples with and without expression data
round_genome_plot_expression <- function(population, expression = T){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palette
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # get input data
  input <- readr::read_csv(paste0("06_local-ancestry-analysis/",
                                  "02_output/",
                                  dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                                                   population == "BORA" ~ "02_BORA",
                                                   population == "BORG" ~ "03_BORG",
                                                   population == "FULA" ~ "04_FULA",
                                                   population == "LAGU" ~ "05_LAGU",
                                                   population == "NDAM" ~ "06_NDAM"),
                                  "/elai-ld-",
                                  population,
                                  dplyr::case_when(expression == T ~ "-expression",
                                                   expression == F ~ ""),
                                  ".csv"))
  # get chromosome lengths
  length <- input %>%
    dplyr::count(chr)
  # set label positions
  label_pos <- c()
  for(i in 1:29){
    label_pos <- c(label_pos,
                   sum(length$n[0:(i-1)],
                       length$n[i]/2))
  }
  # plot results
  plot <- ggplot2::ggplot(input,
                          ggplot2::aes(x = row_num)) +
    ggplot2::geom_area(ggplot2::aes(y = mean_a +
                                      mean_b +
                                      mean_c,
                                    fill = "a")) +
    ggplot2::geom_area(ggplot2::aes(y = mean_b +
                                      mean_c,
                                    fill = "b")) +
    ggplot2::geom_area(ggplot2::aes(y = mean_c,
                                    fill = "c")) +
    ggplot2::annotate("text",
                      x = label_pos,
                      y = -0.1,
                      label = c(1:29),
                      size = c(rep(3,
                                   20),
                               rep(2,
                                   9)),
                      colour = "grey30") +
    ggplot2::scale_fill_manual(labels = c(expression(paste("European ",
                                                           italic("Bos taurus"),
                                                           phantom("p"))),
                                          expression(paste("African ",
                                                           italic("Bos taurus"),
                                                           phantom("p"))),
                                          expression(paste(italic("Bos indicus"),
                                                           phantom("p")))),
                               values = c(pal[17],
                                          pal[9],
                                          pal[28])) +
    ggplot2::scale_x_continuous(expand=c(0,
                                         0)) +
    ggplot2::scale_y_continuous(expand=c(0,
                                         0),
                                limits = c(-1,
                                           1.01)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.text.align = 0,
                   legend.title = ggplot2::element_blank()) +
    ggplot2::coord_polar()
  # add line segments between chromosomes
  for(i in 1:29){
    plot <- plot +
      ggplot2::geom_segment(x = sum(length$n[1:i],
                                    1),
                            xend = sum(length$n[1:i],
                                       1),
                            y = 0,
                            yend = 1,
                            colour = "grey70",
                            linewidth = ggplot2::rel(0.25))
  }
  # save plot
  png(paste0("06_local-ancestry-analysis/",
             "03_figures/",
             dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                              population == "BORA" ~ "02_BORA",
                              population == "BORG" ~ "03_BORG",
                              population == "FULA" ~ "04_FULA",
                              population == "LAGU" ~ "05_LAGU",
                              population == "NDAM" ~ "06_NDAM"),
             "/elai-ld-",
             population,
             dplyr::case_when(expression == T ~ "-expression",
                              expression == F ~ ""),
             "-round.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("06_local-ancestry-analysis/",
             "03_figures/",
             dplyr::case_when(population == "BAOU" ~ "01_BAOU",
                              population == "BORA" ~ "02_BORA",
                              population == "BORG" ~ "03_BORG",
                              population == "FULA" ~ "04_FULA",
                              population == "LAGU" ~ "05_LAGU",
                              population == "NDAM" ~ "06_NDAM"),
             "/elai-ld-",
             population,
             dplyr::case_when(expression == T ~ "-expression",
                              expression == F ~ ""),
             "-round.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

### apply functions to data
#### set populations
populations <- c("BAOU",
                 "BORA",
                 "BORG",
                 "FULA",
                 "LAGU",
                 "NDAM")

#### chromosome plots for all populations with and without expression data
for (i in populations){
  parallel::mclapply(1:29,
                     chromosome_plot_expression,
                     population = i,
                     mc.cores = 15)
  parallel::mclapply(1:29,
                     chromosome_plot_expression,
                     population = i,
                     expression = F,
                     mc.cores = 15)
}

#### flat genome plots for all all populations with and without expression data
parallel::mclapply(populations,
                   flat_genome_plot_expression,
                   mc.cores = 15)
parallel::mclapply(populations,
                   flat_genome_plot_expression,
                   expression = F,
                   mc.cores = 15)

#### round genome plots for all all populations with and without expression data
parallel::mclapply(populations,
                   round_genome_plot_expression,
                   mc.cores = 15)
parallel::mclapply(populations,
                   round_genome_plot_expression,
                   expression = F,
                   mc.cores = 15)