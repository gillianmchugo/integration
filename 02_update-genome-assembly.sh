# update genome assembly
## update genome assembly of snp data to ars-ucd1.2
### remove "-" from snps without identified reference allele in downloaded schnabel file with sed
sed 's/\t-/\t/g' 01_data-sources/06_schnabel/UMC_marker_names_180910/9913_ARS1.2_58336_SNP50_marker_name_180910.map > 01_data-sources/06_schnabel/9913_ARS1.2_58336_SNP50_marker_name_180910-edit.map

### update genome assembly using edited file with plink
#### berthier
plink --cow --bfile 01_data-sources/03_berthier/berthier --update-chr 01_data-sources/06_schnabel/9913_ARS1.2_58336_SNP50_marker_name_180910-edit.map 1 2 --update-map 01_data-sources/06_schnabel/9913_ARS1.2_58336_SNP50_marker_name_180910-edit.map 4 2 --make-bed --out 02_update-genome-assembly/berthier-ars

#### widde
plink --cow --bfile 01_data-sources/04_widde/widde --update-chr 01_data-sources/06_schnabel/9913_ARS1.2_58336_SNP50_marker_name_180910-edit.map 1 2 --update-map 01_data-sources/06_schnabel/9913_ARS1.2_58336_SNP50_marker_name_180910-edit.map 4 2 --make-bed --out 02_update-genome-assembly/widde-ars