# rna-seq analysis
## adapted from AATTOL_R_script_RNAseq_STAR108_EgdeR.r from https://doi.org/10.18167/DVN1/L9SHAX
### install required packages
if (!require(BiocManager)) install.packages("BiocManager")
if (!require(devtools)) install.packages("devtools")
if (!require(edgeR)) BiocManager::install("edgeR")
if (!require(ggpattern)) install.packages("ggpattern")

### import data
#### load count matrix
tab1 <- read.table("01_data-sources/12_gene-expression-omnibus/GSE197108_AATTOL_counts_STAR_120_unik.txt",
                   sep = "\t",
                   header = T,
                   stringsAsFactors = T)

#### load sampling design
samples <- read.table("01_data-sources/11_peylhard/dataverse_files/AATTOL_Samples1_120_star.csv",
                      sep = ";",
                      header = T,
                      stringsAsFactors = T)

#### check dimensions, should be 24616 126
dim(tab1)

#### check dimensions, should be 120 10
dim(samples)

### extract gene counts
#### extract gene info
gene <- tab1[, 1:6]

#### extract gene counts
bos <- tab1[, -c(1:6)]

#### add gene ids
rownames(bos) <- gene$Geneid

### update population names
#### add new column
samples$Breed3 <- samples$Breed

#### set names
levels(samples$Breed3)[1] <- "LAGU"
levels(samples$Breed3)[2] <- "BAOU"
levels(samples$Breed3)[3] <- "BORG"
levels(samples$Breed3)[4] <- "NDAM"
levels(samples$Breed3)[5] <- "FULA"

### add treatment group
samples$Treat <- factor(paste(samples$Breed3,
                              samples$DPI2,
                              sep = "."))

### remove samples Z4, BO5 (outliers) and BA3 (tested positive for trypanosomes before infection)
#### select samples
cZ4BO5BA3 <- c(samples[which(samples$Ani == "Z4"),]$ordre_linux,
               samples[which(samples$Ani == "BO5"),]$ordre_linux,
               samples[which(samples$Ani == "BA3"),]$ordre_linux)

#### remove from sample information
samp_108 <- samples[-cZ4BO5BA3,]

#### remove from data
bos_108 <- bos[,-cZ4BO5BA3]

### create DGElist object without selected samples
DGE_list_bos_108_0 <- edgeR::DGEList(counts = bos_108,
                                     group = samp_108$Treat)

### filter lowly expressed genes
#### count number of samples with count per million (cpm) > 1
A <- rowSums(edgeR::cpm(DGE_list_bos_108_0) > 1)

#### keep libraries with at least two samples with cpm > 1
keep <- rowSums(edgeR::cpm(DGE_list_bos_108_0) > 1) >= 2

#### check summary, should be 11509 false 13107 true
summary(keep)

#### apply filter to data
DGEL_bos_108_F <- DGE_list_bos_108_0[keep, ,
                                     keep.lib.sizes = F]

#### check dimensions, should be 13107 108
dim(DGEL_bos_108_F$counts)

#### get total number of genes
geneF <- dim(DGEL_bos_108_F$counts)[1]

### normalize library sizes
DGEL_bos_108_F <- edgeR::calcNormFactors(DGEL_bos_108_F)

### de analysis
#### create design matrix
response_design <- model.matrix(~ 0 + samp_108$Treat)

#### extract and modify column names
response_a <- substr(colnames(response_design),
                     10,
                     21)

#### replace column names
colnames(response_design) <- response_a

#### create contrasts
response_contrasts <- limma::makeContrasts(LAGU_20 = (TreatLAGU.20 - TreatLAGU.0) - (TreatFULA.20 - TreatFULA.0),
                                           LAGU_30 = (TreatLAGU.30 - TreatLAGU.0) - (TreatFULA.30 - TreatFULA.0),
                                           LAGU_40 = (TreatLAGU.40 - TreatLAGU.0) - (TreatFULA.40 - TreatFULA.0),
                                           NDAM_20 = (TreatNDAM.20 - TreatNDAM.0) - (TreatFULA.20 - TreatFULA.0),
                                           NDAM_30 = (TreatNDAM.30 - TreatNDAM.0) - (TreatFULA.30 - TreatFULA.0),
                                           NDAM_40 = (TreatNDAM.40 - TreatNDAM.0) - (TreatFULA.40 - TreatFULA.0),
                                           BAOU_20 = (TreatBAOU.20 - TreatBAOU.0) - (TreatFULA.20 - TreatFULA.0),
                                           BAOU_30 = (TreatBAOU.30 - TreatBAOU.0) - (TreatFULA.30 - TreatFULA.0),
                                           BAOU_40 = (TreatBAOU.40 - TreatBAOU.0) - (TreatFULA.40 - TreatFULA.0),
                                           BORG_20 = (TreatBORG.20 - TreatBORG.0) - (TreatFULA.20 - TreatFULA.0),
                                           BORG_30 = (TreatBORG.30 - TreatBORG.0) - (TreatFULA.30 - TreatFULA.0),
                                           BORG_40 = (TreatBORG.40 - TreatBORG.0) - (TreatFULA.40 - TreatFULA.0),
                                           levels = response_design)

#### estimate dispersions
response_DGEL_bos_108_F <- edgeR::estimateDisp(DGEL_bos_108_F,
                                               response_design)

#### fit general linear model
response_DGE_fit_bos_108 <- edgeR::glmFit(response_DGEL_bos_108_F,
                                          response_design)

#### function to extract results
f_test_edger2 <- function(fit, contraste){
  # genewise likelihood ratio tests
  tmp <- edgeR::glmLRT(fit,
                       contrast = contraste)
  # top de genes
  tmp2 <- edgeR::topTags(tmp,
                         n = geneF)
  tmp3 <- tmp2$table
  res <- tmp3[order(rownames(tmp3)),]
}

#### apply function to extract results
toptags_LAGU_20 <- f_test_edger2(response_DGE_fit_bos_108,
                                 response_contrasts[,"LAGU_20"])
toptags_LAGU_30 <- f_test_edger2(response_DGE_fit_bos_108,
                                 response_contrasts[,"LAGU_30"])
toptags_LAGU_40 <- f_test_edger2(response_DGE_fit_bos_108,
                                 response_contrasts[,"LAGU_40"])
toptags_NDAM_20 <- f_test_edger2(response_DGE_fit_bos_108,
                                 response_contrasts[,"NDAM_20"])
toptags_NDAM_30 <- f_test_edger2(response_DGE_fit_bos_108,
                                 response_contrasts[,"NDAM_30"])
toptags_NDAM_40 <- f_test_edger2(response_DGE_fit_bos_108,
                                 response_contrasts[,"NDAM_40"])
toptags_BAOU_20 <- f_test_edger2(response_DGE_fit_bos_108,
                                 response_contrasts[,"BAOU_20"])
toptags_BAOU_30 <- f_test_edger2(response_DGE_fit_bos_108,
                                 response_contrasts[,"BAOU_30"])
toptags_BAOU_40 <- f_test_edger2(response_DGE_fit_bos_108,
                                 response_contrasts[,"BAOU_40"])
toptags_BORG_20 <- f_test_edger2(response_DGE_fit_bos_108,
                                 response_contrasts[,"BORG_20"])
toptags_BORG_30 <- f_test_edger2(response_DGE_fit_bos_108,
                                 response_contrasts[,"BORG_30"])
toptags_BORG_40 <- f_test_edger2(response_DGE_fit_bos_108,
                                 response_contrasts[,"BORG_40"])

#### combine results
response_edger_13107 <- cbind(toptags_LAGU_20[, c(1,5)],
                              toptags_LAGU_30[, c(1, 5)],
                              toptags_LAGU_40[, c(1, 5)],
                              toptags_NDAM_20[, c(1, 5)],
                              toptags_NDAM_30[, c(1, 5)],
                              toptags_NDAM_40[, c(1, 5)],
                              toptags_BAOU_20[, c(1, 5)],
                              toptags_BAOU_30[, c(1, 5)],
                              toptags_BAOU_40[, c(1, 5)],
                              toptags_BORG_20[, c(1, 5)],
                              toptags_BORG_30[, c(1, 5)],
                              toptags_BORG_40[, c(1, 5)])

#### add id column
response_edger_13107$ID <- row.names(response_edger_13107)

#### add column names
colnames(response_edger_13107) <- c("logFC_LAGU_20", "FDR_LAGU_20",
                                    "logFC_LAGU_30", "FDR_LAGU_30",
                                    "logFC_LAGU_40", "FDR_LAGU_40",
                                    "logFC_NDAM_20", "FDR_NDAM_20",
                                    "logFC_NDAM_30", "FDR_NDAM_30",
                                    "logFC_NDAM_40", "FDR_NDAM_40",
                                    "logFC_BAOU_20", "FDR_BAOU_20",
                                    "logFC_BAOU_30", "FDR_BAOU_30",
                                    "logFC_BAOU_40", "FDR_BAOU_40",
                                    "logFC_BORG_20", "FDR_BORG_20",
                                    "logFC_BORG_30", "FDR_BORG_30",
                                    "logFC_BORG_40", "FDR_BORG_40",
                                    "ID")

#### load annotation file
gene_annot_13107_E92 <- read.table("01_data-sources/10_peylhard/dataverse_files/gene_annot_13107_IPA_E92.1.txt",
                                   sep = "\t",
                                   header = T)

#### check summary, should all be true
summary(as.character(rownames(response_edger_13107)) == as.character(gene_annot_13107_E92$ID))

#### merge results and annotation with ensembl ids
tmp <- response_edger_13107
tmp$ID <- rownames(tmp)
tmp <- merge(tmp,
             gene_annot_13107_E92[, c(1, 18)])
response_edger_13107_IN_IPA <- tmp

#### write results with ensembl ids
write.table(response_edger_13107_IN_IPA,
            file = "07_differential-expression-analysis/01_output/response_edger_13107_IN_IPA.txt",
            quote = F,
            row.names = F,
            sep = "\t")

#### merge, save and read results and annotation with gene symbols
tmp <- response_edger_13107
tmp$ID <- rownames(tmp)
tmp <- merge(tmp,
             gene_annot_13107_E92[, c(1, 2)])
response_rna_seq_gene_symbols <- tmp
readr::write_csv(response_rna_seq_gene_symbols,
                 file = "07_differential-expression-analysis/01_output/response-rna-seq-gene-symbols.csv")
response_rna_seq_gene_symbols <- readr::read_csv("07_differential-expression-analysis/01_output/response-rna-seq-gene-symbols.csv")

#### set contrast ids
rna_seq_response_contrast_ids <- c("BORG_20",
                                   "NDAM_20",
                                   "BAOU_20",
                                   "LAGU_20",
                                   "BORG_30",
                                   "NDAM_30",
                                   "BAOU_30",
                                   "LAGU_30",
                                   "BORG_40",
                                   "NDAM_40",
                                   "BAOU_40",
                                   "LAGU_40")

#### extract all de symbols
all_response_rna_seq_de_gene_symbols <- data.frame(Symbol = NA,
                                                   set = NA)
for (i in 1:length(rna_seq_response_contrast_ids)){
  all_increased_response_rna_seq_de_gene_symbols <- response_rna_seq_gene_symbols %>%
    dplyr::filter(!!rlang::sym(paste0("FDR_",
                                      rna_seq_response_contrast_ids[i])) <=  0.05) %>%
    dplyr::filter(!!rlang::sym(paste0("logFC_",
                                      rna_seq_response_contrast_ids[i])) >  0) %>%
    dplyr::arrange(!!rlang::sym(paste0("FDR_",
                                       rna_seq_response_contrast_ids[i]))) %>%
    dplyr::distinct(Symbol,
                    .keep_all = T) %>%
    dplyr::filter(Symbol != " ") %>%
    dplyr::select(Symbol) %>%
    dplyr::mutate(set = paste0(rna_seq_response_contrast_ids[i],
                               "_increased"))
  all_decreased_response_rna_seq_de_gene_symbols <- response_rna_seq_gene_symbols %>%
    dplyr::filter(!!rlang::sym(paste0("FDR_",
                                      rna_seq_response_contrast_ids[i])) <=  0.05) %>%
    dplyr::filter(!!rlang::sym(paste0("logFC_",
                                      rna_seq_response_contrast_ids[i])) <  0) %>%
    dplyr::arrange(!!rlang::sym(paste0("FDR_",
                                       rna_seq_response_contrast_ids[i]))) %>%
    dplyr::distinct(Symbol,
                    .keep_all = T) %>%
    dplyr::filter(Symbol != " ") %>%
    dplyr::select(Symbol) %>%
    dplyr::mutate(set = paste0(rna_seq_response_contrast_ids[i],
                               "_decreased"))
  all_response_rna_seq_de_gene_symbols <- dplyr::full_join(all_response_rna_seq_de_gene_symbols,
                                                           all_increased_response_rna_seq_de_gene_symbols)
  all_response_rna_seq_de_gene_symbols <- dplyr::full_join(all_response_rna_seq_de_gene_symbols,
                                                           all_decreased_response_rna_seq_de_gene_symbols)
}
response_rna_seq_de_numbers <- all_response_rna_seq_de_gene_symbols %>%
  tidyr::drop_na() %>%
  dplyr::group_by(set) %>%
  dplyr::summarise(number = dplyr::n()) %>%
  dplyr::mutate(expression = stringr::str_split_i(set,
                                                  "_",
                                                  3),
                dpi = stringr::str_split_i(set,
                                           "_",
                                           2),
                population = stringr::str_split_i(set,
                                                  "_",
                                                  1),
                contrast_id = paste0(population,
                                     "_",
                                     dpi)) %>%
  dplyr::arrange(dpi,
                 population,
                 desc(expression)) %>%
  readr::write_csv("07_differential-expression-analysis/01_output/rna-seq-de-gene-numbers.csv")

response_rna_seq_de_numbers <- readr::read_csv("07_differential-expression-analysis/01_output/rna-seq-de-gene-numbers.csv")

### function to generate colour palette
smooth_rainbow <- khroma::colour("smooth rainbow")

###generate colour palettes
pal <- smooth_rainbow(31,
                      range = c(0.1,
                                0.9))

pal2 <- khroma::colour("sunset")(11)

viridis <- viridis::viridis(41, direction = -1)

### set colour palettes
cols <-  c(pal[8],
           pal2[[1]],
           pal[5],
           pal[4])

dpi_cols <- c(viridis[21],
              viridis[31],
              viridis[41])

### set limits
lims <- c("LAGU",
          "BAOU",
          "NDAM",
          "BORG")

### draw bar plot
plot <- ggplot2::ggplot(response_rna_seq_de_numbers,
                        ggplot2::aes(x = factor(interaction(population,
                                                            dpi),
                                                levels = c("LAGU.20",
                                                           "BAOU.20",
                                                           "NDAM.20",
                                                           "BORG.20",
                                                           "LAGU.30",
                                                           "BAOU.30",
                                                           "NDAM.30",
                                                           "BORG.30",
                                                           "LAGU.40",
                                                           "BAOU.40",
                                                           "NDAM.40",
                                                           "BORG.40")),
                                     colour = population,
                                     fill = population,
                                     pattern_shape = population)) +
  ggpattern::geom_col_pattern(ggplot2::aes(y = number),
                              data = subset(response_rna_seq_de_numbers,
                                            expression == "increased"),
                              linewidth = 1,
                              pattern = "pch",
                              pattern_alpha = 0.5,
                              pattern_angle = 0,
                              pattern_colour = NA,
                              pattern_density = 0.6,
                              pattern_fill = "white",
                              pattern_spacing = 0.03) +
  ggpattern::geom_col_pattern(ggplot2::aes(y = -number),
                              data = subset(response_rna_seq_de_numbers,
                                            expression == "decreased"),
                              linewidth = 1,
                              pattern = "pch",
                              pattern_alpha = 0.5,
                              pattern_angle = 0,
                              pattern_colour = NA,
                              pattern_density = 0.6,
                              pattern_fill = "white",
                              pattern_spacing = 0.03) +
  ggplot2::annotate("text",
                    label = c("Increased Expression",
                              "Decreased Expression"),
                    hjust = -0.05,
                    x = -Inf,
                    y = c(21*0.917,
                          -21*0.917)) +
  ggplot2::geom_hline(colour = "grey70",
                      linewidth = ggplot2::rel(0.25),
                      yintercept = 0, ) +
  ggplot2::scale_colour_manual(name = "Population",
                               limits = lims,
                               values = cols) +
  ggplot2::scale_fill_manual(name = "Population",
                             limits = lims,
                             values = cols) +
  ggpattern::scale_pattern_shape_manual(limits = lims,
                                        name = "Population",
                                        values = c(22,
                                                   24,
                                                   21,
                                                   22)) +
  ggplot2::scale_x_discrete(drop = F,
                            guide = ggh4x::guide_axis_nested()) +
  ggplot2::scale_y_continuous(expand = c(0,
                                         0),
                              limits = c(-21,
                                         21)) +
  ggplot2::labs(x = "Days Post-Infection",
                y = "Number of Significantly Differentially Expressed Genes") +
  ggplot2::theme_light() +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 ggh4x.axis.nesttext.x = ggplot2::element_text(),
                 ggh4x.axis.nestline = ggplot2::element_line(colour = dpi_cols,
                                                             linewidth = 0.5),
                 legend.position = "top",
                 panel.border = ggplot2::element_rect(linewidth = ggplot2::rel(0.5)))

### save bar plot
png("07_differential-expression-analysis/02_figures/bar-plot.png",
    height = 6,
    res = 300,
    units = "in",
    width = 9)
print(plot)
dev.off()
pdf("07_differential-expression-analysis/02_figures/bar-plot.pdf",
    height = 6,
    width = 9)
print(plot)
dev.off()

### save all gene symbols
all_response_rna_seq_de_gene_symbols %>%
  tidyr::drop_na() %>%
  dplyr::group_by(set) %>%
  dplyr::summarise(genes = toString(Symbol)) %>%
  dplyr::mutate(expression = stringr::str_split_i(set,
                                                  "_",
                                                  3),
                dpi = stringr::str_split_i(set,
                                           "_",
                                           2),
                population = stringr::str_split_i(set,
                                                  "_",
                                                  1)) %>%
  dplyr::arrange(population,
                 desc(expression),
                 dpi) %>%
  dplyr::relocate(genes,
                  .after = last_col()) %>%
  readr::write_csv("07_differential-expression-analysis/01_output/all-rna-seq-de-gene-symbols.csv")

### extract top de gene symbols
top_response_rna_seq_de_gene_symbols <- data.frame(Symbol = NA, set = NA)
for (i in 1:length(rna_seq_response_contrast_ids)){
  top_increased_response_rna_seq_de_gene_symbols <- response_rna_seq_gene_symbols %>%
    dplyr::filter(!!rlang::sym(paste0("FDR_",
                                      rna_seq_response_contrast_ids[i])) <=  0.05) %>%
    dplyr::filter(!!rlang::sym(paste0("logFC_",
                                      rna_seq_response_contrast_ids[i])) >  0) %>%
    dplyr::distinct(Symbol,
                    .keep_all = T) %>%
    dplyr::filter(Symbol != " ") %>%
    dplyr::slice_min(n = 10,
                     order_by = !!rlang::sym(paste0("FDR_",
                                                    rna_seq_response_contrast_ids[i]))) %>%
    dplyr::select(Symbol) %>%
    dplyr::mutate(set = paste0(rna_seq_response_contrast_ids[i],
                               "_increased"))
  top_decreased_response_rna_seq_de_gene_symbols <- response_rna_seq_gene_symbols %>%
    dplyr::filter(!!rlang::sym(paste0("FDR_",
                                      rna_seq_response_contrast_ids[i])) <=  0.05) %>%
    dplyr::filter(!!rlang::sym(paste0("logFC_",
                                      rna_seq_response_contrast_ids[i])) < 0) %>%
    dplyr::distinct(Symbol,
                    .keep_all = T) %>%
    dplyr::filter(Symbol != " ") %>%
    dplyr::slice_min(n = 10,
                     order_by = !!rlang::sym(paste0("FDR_",
                                                    rna_seq_response_contrast_ids[i]))) %>%
    dplyr::select(Symbol) %>%
    dplyr::mutate(set = paste0(rna_seq_response_contrast_ids[i],
                               "_decreased"))
  top_response_rna_seq_de_gene_symbols <- dplyr::full_join(top_response_rna_seq_de_gene_symbols,
                                                           top_increased_response_rna_seq_de_gene_symbols)
  top_response_rna_seq_de_gene_symbols <- dplyr::full_join(top_response_rna_seq_de_gene_symbols,
                                                           top_decreased_response_rna_seq_de_gene_symbols)
}
top_response_rna_seq_de_gene_symbols %>%
  tidyr::drop_na() %>%
  dplyr::group_by(set) %>%
  dplyr::summarise(genes = toString(Symbol)) %>%
  dplyr::mutate(expression = stringr::str_split_i(set,
                                                  "_",
                                                  3),
                dpi = stringr::str_split_i(set,
                                           "_",
                                           2),
                population = stringr::str_split_i(set,
                                                  "_",
                                                  1)) %>%
  dplyr::arrange(population,
                 desc(expression),
                 dpi) %>%
  dplyr::relocate(genes,
                  .after = last_col()) %>%
  readr::write_csv("07_differential-expression-analysis/01_output/top-response-rna-seq-de-gene-symbols.csv")

### function to draw and save volcano plots
volcano_plot <- function(contrast){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palettes
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # set colour palettes
  cols <-  c(pal[26],
             "#666666",
             pal[5])
  # set labels
  labs <- c("Significantly Decreased Expression",
            "No Significant Change",
            "Significantly Increased Expression")
  # set limits
  lims <- c(-1,
            0,
            1)
  # prepare input
  input <- response_rna_seq_gene_symbols %>%
    dplyr::mutate(results = dplyr::case_when(!!rlang::sym(paste0("FDR_",
                                                                 contrast)) <= 0.05 &
                                               !!rlang::sym(paste0("logFC_",
                                                                   contrast)) > 0 ~ 1,
                                             !!rlang::sym(paste0("FDR_",
                                                                 contrast)) <= 0.05 &
                                               !!rlang::sym(paste0("logFC_",
                                                                   contrast)) < 0 ~ -1,
                                             .default = 0))
  # select top 10 genes
  top_10_up <- input %>%
    tidyr::drop_na(Symbol) %>%
    dplyr::filter(Symbol != " ") %>%
    dplyr::filter(results == 1) %>%
    dplyr::arrange(!!rlang::sym(paste0("FDR_",
                                       contrast))) %>%
    dplyr::distinct(Symbol,
                    .keep_all = T) %>%
    dplyr::slice_min(n = 10,
                     order_by = !!rlang::sym(paste0("FDR_",
                                                    contrast)))
  top_10_down <- input %>%
    tidyr::drop_na(Symbol) %>%
    dplyr::filter(Symbol != " ") %>%
    dplyr::filter(results == -1) %>%
    dplyr::arrange(!!rlang::sym(paste0("FDR_",
                                       contrast))) %>%
    dplyr::distinct(Symbol,
                    .keep_all = T) %>%
    dplyr::slice_min(n = 10,
                     order_by = !!rlang::sym(paste0("FDR_",
                                                    contrast)))
  top_20 <- dplyr::full_join(top_10_up,
                             top_10_down)
  # set x max
  x_max <- input %>%
    dplyr::slice_max(n = 1,
                     with_ties = F,
                     order_by = abs(!!rlang::sym(paste0("logFC_",
                                                        contrast)))) %>%
    dplyr::pull(!!rlang::sym(paste0("logFC_",
                                    contrast))) %>%
    abs()
  # draw plot
  plot <- input %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[paste0("logFC_",
                                                   contrast)]],
                                 y = -log(.data[[paste0("FDR_",
                                                        contrast)]],
                                          base = 10),
                                 colour = as.factor(results),
                                 fill = as.factor(results))) +
    ggplot2::geom_point(shape = 21) +
    ggplot2::geom_hline(colour = "grey70",
                        linetype = "dashed",
                        yintercept = -log(0.05,
                                          base = 10)) +
    ggplot2::labs(x = "Log<sub>2</sub> Fold Change",
                  y = "-Log<sub>10</sub>*P*<sub>adj.</sub>") +
    ggplot2::scale_colour_manual(labels = labs,
                                 limits = as.factor(lims),
                                 values = ggplot2::alpha(cols, c(1,
                                                                 0.1,
                                                                 1))) +
    ggplot2::scale_fill_manual(labels = labs,
                               limits = as.factor(lims),
                               values = ggplot2::alpha(cols, c(0.5,
                                                               0.05,
                                                               0.5))) +
    ggplot2::scale_x_continuous(limits = c(-x_max,
                                           x_max)) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.title.x = ggtext::element_markdown(),
                   axis.title.y = ggtext::element_markdown(),
                   legend.position = "top",
                   legend.title = ggplot2::element_blank())
  # save plot
  png(paste0("07_differential-expression-analysis/02_figures/",
             stringr::str_replace_all(contrast,
                                      "_",
                                      "-"),
             "-volcano.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("07_differential-expression-analysis/02_figures/",
             stringr::str_replace_all(contrast,
                                      "_",
                                      "-"),
             "-volcano.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
  # add labels
  plot_labels <- plot +
    ggrepel::geom_text_repel(data = top_20,
                             ggplot2::aes(label = Symbol),
                             bg.color = "white",
                             bg.r = 0.1,
                             box.padding = 0.5,
                             fontface = "italic",
                             force = 20,
                             force_pull = 0,
                             lineheight = 0.7,
                             max.iter = 1000000000,
                             max.overlaps = Inf,
                             max.time = 10,
                             min.segment.length = 0,
                             segment.size = 0.25,
                             show.legend = F,
                             size = ggplot2::rel(3))
  # save plot
  png(paste0("07_differential-expression-analysis/02_figures/",
             stringr::str_replace_all(contrast,
                                      "_",
                                      "-"),
             "-volcano-labels.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot_labels)
  dev.off()
  pdf(paste0("07_differential-expression-analysis/02_figures/",
             stringr::str_replace_all(contrast,
                                      "_",
                                      "-"),
             "-volcano-labels.pdf"),
      height = 6,
      width = 9)
  print(plot_labels)
  dev.off()
}

### apply volcano plot function
lapply(rna_seq_response_contrast_ids,
       volcano_plot)