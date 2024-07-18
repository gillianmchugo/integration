# structure analysis
## structure analysis with faststructure
### run structure analysis up to k = 25 with faststructure using structure_threader
nohup structure_threader run -fs /home/workspace/gmchugo/software/fastStructure/structure.py -t 5 -K 27 --no_plots 1 -i 03_filter-data/05_filter-snps/ld.bed --pop 05_structure-analysis/01_input/popfile.txt -o 05_structure-analysis/02_output > 05_structure-analysis/02_output/log.txt &

### plot results
Rscript 05_structure-analysis/03_figures/structure-plots.R
Rscript 05_structure-analysis/03_figures/structure-plots-labels.R