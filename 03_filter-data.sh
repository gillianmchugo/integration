# filter data
## filter snp data using plink and r
### merge files with plink
plink --cow --bfile 02_update-genome-assembly/berthier-ars --bmerge 02_update-genome-assembly/widde-ars --make-bed --out 03_filter-data/01_low-density-snps/berthier-widde

### extract low density snps
plink --cow --bfile 03_filter-data/01_low-density-snps/berthier-widde --recode --out 03_filter-data/01_low-density-snps/berthier-widde
Rscript 03_filter-data/01_low-density-snps/ld-snps.R
plink --cow --bfile 01_data-sources/01_unfiltered-snp-data/ars-mind-ibs-het-filter --extract 03_filter-data/01_low-density-snps/ld-snps.txt --make-bed --out 03_filter-data/01_low-density-snps/ld-snps

### merge files with plink
plink --cow --bfile 03_filter-data/01_low-density-snps/berthier-widde --bmerge 03_filter-data/01_low-density-snps/ld-snps --make-bed --out 03_filter-data/01_low-density-snps/merged

### missing genotype filter to remove individuals with missing ld snps
plink --cow --bfile 03_filter-data/01_low-density-snps/merged --mind 0.95 --make-bed --out 03_filter-data/02_missing-snps/merged

### identity by state filter to remove duplicate individuals
plink --cow --bfile 03_filter-data/01_low-density-snps/merged --distance ibs square --out 03_filter-data/03_identity-by-state/merged
Rscript 03_filter-data/03_identity-by-state/ibs-plot.R
Rscript 03_filter-data/03_identity-by-state/ibs-filter.R
plink --cow --bfile 03_filter-data/01_low-density-snps/merged --remove 03_filter-data/03_identity-by-state/ibs-filter.txt --make-bed --out 03_filter-data/03_identity-by-state/merged-ibs

### inbreeding
plink --cow --bfile 03_filter-data/03_identity-by-state/merged-ibs --het --out 03_filter-data/04_inbreeding/merged-ibs
Rscript 03_filter-data/04_inbreeding/inbreeding-plot.R

### filter snps by autosome, call rate and minor allele frequency
plink --cow --bfile 03_filter-data/03_identity-by-state/merged-ibs --chr 1-29 --geno 0.05 --maf 0.05 --make-bed --out 03_filter-data/05_filter-snps/ld
