# local ancestry analysis
## local ancestry analysis with elai
### generate bimbam format files for each breed and chromosome using plink
nohup bash -c 'for b in ALEN ANGU ANKO BAOU BORA BORG CHIA EASZ FULA GIR HOLS JERS KARA KETE LAGU MARC MARE MUTU NDAG NDAM NELO NGAN ROMA SHEK SOMB THAR; do for i in {1..29}; do plink --cow --bfile 03_filter-data/05_filter-snps/ld --family --keep-cluster-names $b --chr $i --recode-bimbam --out 06_local-ancestry-analysis/01_input/ld-$b-$i; done; done' > 06_local-ancestry-analysis/01_input/log.txt &

### run elai analysis
cd 06_local-ancestry-analysis/02_output/01_BAOU
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-NDAG-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-BAOU-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-BAOU-$i.recode.pos.txt -s 30 -o ld-BAOU-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 02_BORA
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-NDAG-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-BORA-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-BORA-$i.recode.pos.txt -s 30 -o ld-BORA-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 03_BORG
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-NDAG-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-BORG-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-BORG-$i.recode.pos.txt -s 30 -o ld-BORG-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 04_FULA
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-NDAG-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-FULA-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-FULA-$i.recode.pos.txt -s 30 -o ld-FULA-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 05_LAGU
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-NDAG-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-LAGU-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-LAGU-$i.recode.pos.txt -s 30 -o ld-LAGU-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 06_NDAM
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-NDAG-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-NDAM-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/integration/06_local-ancestry-analysis/01_input/ld-NDAM-$i.recode.pos.txt -s 30 -o ld-NDAM-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ../../..

### extract local ancestry results
nohup bash -c '06_local-ancestry-analysis/02_output/elai-results.R' > /dev/null 2>&1 &

### generate local ancestry plots
nohup bash -c 'Rscript 06_local-ancestry-analysis/03_figures/local-ancestry-plots.R' > 06_local-ancestry-analysis/03_figures/log.txt &

### combine local ancestry plots
nohup bash -c 'Rscript 06_local-ancestry-analysis/03_figures/combine-local-ancestry-plots.R' > /dev/null 2>&1 &
