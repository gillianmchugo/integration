# integration analysis
## integration analysis using inrich
### prepare integration input
nohup bash -c '10_integration-analysis/01_input/integration-input.R' > /dev/null 2>&1 &

### generate clumped files with plink
nohup bash -c 'for i in 10_integration-analysis/01_input/04_intervals/*hd*hybrids-plink-input.txt; do for j in a b c; do plink --cow --bfile 01_data-sources/02_filtered-snp-data/hd --clump $i --clump-snp-field rs --clump-field p_$j --clump-p1 0.02275 --clump-p2 0.05 --clump-kb 1000 --out $i-$j; done; done' > 10_integration-analysis/01_input/04_intervals/log.txt &
nohup bash -c 'for i in 10_integration-analysis/01_input/04_intervals/*ld*hybrids-plink-input.txt; do for j in a b c; do plink --cow --bfile 01_data-sources/02_filtered-snp-data/ld --clump $i --clump-snp-field rs --clump-field p_$j --clump-p1 0.02275 --clump-p2 0.05 --clump-kb 1000 --out $i-$j; done; done' > 10_integration-analysis/01_input/04_intervals/log.txt &
nohup bash -c 'for i in 10_integration-analysis/01_input/04_intervals/*hd*microarray-plink-input.txt; do for j in a b c; do plink --cow --bfile 01_data-sources/02_filtered-snp-data/hd --clump $i --clump-snp-field rs --clump-field p_$j --clump-p1 0.02275 --clump-p2 0.05 --clump-kb 1000 --out $i-$j; done; done' > 10_integration-analysis/01_input/04_intervals/log.txt &
nohup bash -c 'for i in 10_integration-analysis/01_input/04_intervals/*ld*microarray-plink-input.txt; do for j in a b c; do plink --cow --bfile 01_data-sources/02_filtered-snp-data/ld --clump $i --clump-snp-field rs --clump-field p_$j --clump-p1 0.02275 --clump-p2 0.05 --clump-kb 1000 --out $i-$j; done; done' > 10_integration-analysis/01_input/04_intervals/log.txt &
nohup bash -c 'for i in 10_integration-analysis/01_input/04_intervals/*expression-plink-input.txt; do for j in a b c; do plink --cow --bfile 03_filter-data/05_filter-snps/ld --clump $i --clump-snp-field rs --clump-field p_$j --clump-p1 0.02275 --clump-p2 0.05 --clump-kb 1000 --out $i-$j; done; done' > 10_integration-analysis/01_input/04_intervals/log.txt &
nohup bash -c 'for i in 10_integration-analysis/01_input/04_intervals/*expression-all-plink-input.txt; do for j in a b c; do plink --cow --bfile 03_filter-data/05_filter-snps/ld --clump $i --clump-snp-field rs --clump-field p_$j --clump-p1 0.02275 --clump-p2 0.05 --clump-kb 1000 --out $i-$j; done; done' > 10_integration-analysis/01_input/04_intervals/log.txt &

### prepare interval files
nohup bash -c '10_integration-analysis/01_input/04_intervals/interval-files.R' > /dev/null 2>&1 &

### move directory
cd 10_integration-analysis

## install inrich
wget https://zzz.bwh.harvard.edu/inrich/bin/inrich.v.1.1.tar.gz
tar xvzf inrich.v.1.1.tar.gz

### perform enrichment analysis
nohup bash -c 'for i in 01_input/04_intervals/*expression*intervals.txt; do ./inrich -a $i -m 01_input/02_snps/merged-ld-snps.txt -g 01_input/01_genes/genes.txt -t 01_input/03_targets/end-module-targets.txt -u -j 2500 -o 02_output/$i; done' > log.txt &
nohup bash -c 'for i in 01_input/04_intervals/*hd*intervals.txt; do ./inrich -a $i -m 01_input/02_snps/hd-snps.txt -g 01_input/01_genes/genes.txt -t 01_input/03_targets/end-module-targets.txt -u -j 2500 -o 02_output/$i; done' > log.txt &
nohup bash -c 'for i in 01_input/04_intervals/*ld*intervals.txt; do ./inrich -a $i -m 01_input/02_snps/ld-snps.txt -g 01_input/01_genes/genes.txt -t 01_input/03_targets/end-module-targets.txt -u -j 2500 -o 02_output/$i; done' > log.txt &

### format output of integration analysis
nohup bash -c '10_integration-analysis/02_output/integration-output.R' > /dev/null 2>&1 &