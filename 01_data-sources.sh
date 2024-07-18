# data sources
## copy or download files with wget or web interfaces
### low density snp data
#### copy from local ancestry folder

### filtered snp data
#### copy from local ancestry folder

### berthier
#### download data from https://doi.org/10.18167/DVN1/APTZOC

### snpchimp
#### select required data from https://webserver.ibba.cnr.it/SNPchimp/index.php/download/download-cow-data

### schnabel
wget https://www.animalgenome.org/repository/cattle/UMC_bovine_coordinates/UMC_marker_names_180910.zip -O 01_data-sources/06_schnabel/UMC_marker_names_180910.zip

### principal component results
#### copy from local ancestry folder

### structure results
#### copy from local ancestry folder

### local ancestry results
#### copy from local ancestry folder

### differential expression results
#### copy from microarray folder

### peylhard
#### download data from https://doi.org/10.18167/DVN1/L9SHAX

### gene expression omnibus
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197108/suppl/GSE197108%5FAATTOL%5Fcounts%5FSTAR%5F120%5Funik%2Etxt%2Egz -O 01_data-sources/12_gene-expression-omnibus/GSE197108_AATTOL_counts_STAR_120_unik.txt.gz

## unzip downloaded files with unzip or gunzip
### berthier
unzip 01_data-sources/03_berthier/dataverse_files.zip -d 01_data-sources/03_berthier/dataverse_files

### widde
unzip 01_data-sources/04_widde/cattle_*_PLINK.zip -d 01_data-sources/04_widde

### snpchimp
gunzip 01_data-sources/05_snpchimp/SNPchimp_result_*.csv.gz

### schnabel
unzip 01_data-sources/06_schnabel/UMC_marker_names_180910.zip -d 01_data-sources/06_schnabel

### peylhard
unzip 01_data-sources/11_peylhard/dataverse_files.zip -d 01_data-sources/11_peylhard/dataverse_files

### gene expression omnibus
gunzip 01_data-sources/12_gene-expression-omnibus/GSE197108_AATTOL_counts_STAR_120_unik.txt.gz

## convert snp data to binary plink files in forward format with plink and snpchimp
### berthier
python /path/to/snpchimp/source_codes/iConvert/iConvert.py -p 01_data-sources/03_berthier/convert.param
plink --cow --ped 01_data-sources/03_berthier/dataverse_files/cattle__54229variants__39individuals_updated.ped --map 01_data-sources/03_berthier/dataverse_files/cattle__54229variants.map --make-bed --out 01_data-sources/03_berthier/cattle__54229variants__39individuals_updated

### widde
python /path/to/snpchimp/source_codes/iConvert/iConvert.py -p 01_data-sources/04_widde/convert.param
plink --cow --ped 01_data-sources/04_widde/cattle__50975variants__370individuals_updated.ped --map 01_data-sources/04_widde/cattle__50975variants__370individuals.map --make-bed --out 01_data-sources/04_widde/cattle__50975variants__370individuals_updated

## update snp data sample ids using edited fam files with plink
### berthier
plink --cow --bfile 01_data-sources/03_berthier/cattle__54229variants__39individuals_updated --update-ids 01_data-sources/03_berthier/berthier-ids.txt --make-bed --out 01_data-sources/03_berthier/berthier

### widde
plink --cow --bfile 01_data-sources/04_widde/cattle__50975variants__370individuals_updated --update-ids 01_data-sources/04_widde/widde-ids.txt --make-bed --out 01_data-sources/04_widde/widde