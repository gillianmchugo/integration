# principal component analysis
## principal component analysis of snp data with eigensoft
### copy bed file to directory with cp
cp 03_filter-data/05_filter-snps/ld.bed 04_principal-component-analysis/01_input

### copy bim file renamed as pedsnp file with cp
cp 03_filter-data/05_filter-snps/ld.bim 04_principal-component-analysis/01_input/ld.pedsnp

### copy fam file renamed as pedind file and replace -9 with 1 with sed
sed 's/-9/1/g' 03_filter-data/05_filter-snps/ld.fam > 04_principal-component-analysis/01_input/ld.pedind

### convert files using parameter file with convertf
/path/to/EIG-7.2.1/src/convertf -p 04_principal-component-analysis/01_input/ld-convertf.txt

### export paths to software dependencies
export LD_LIBRARY_PATH=/path/to/openblas/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/gsl/lib

### run principal component analysis using parameter file with smartpca
/path/to/EIG-7.2.1/src/eigensrc/smartpca -p 04_principal-component-analysis/01_input/ld-smartpca.txt

### generate figures
Rscript 04_principal-component-analysis/03_figures/pca-plots.R