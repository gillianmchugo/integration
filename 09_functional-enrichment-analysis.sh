# functional enrichment analysis
## functional enrichment analysis of local ancestry, de and network results
### functional enrichment analysis of local ancestry results
nohup bash -c '09_functional-enrichment-analysis/01_local-ancestry-analysis/local-ancestry-functional-enrichment.R' > /dev/null 2>&1 &

### combine local ancestry enrichment plots
nohup bash -c 'Rscript 09_functional-enrichment-analysis/01_local-ancestry-analysis/03_figures/combine-local-ancestry-enrichment-plots.R' > /dev/null 2>&1 &

### functional enrichment analysis of de results
nohup bash -c '09_functional-enrichment-analysis/02_differential-expression-analysis/differential-expression-functional-enrichment.R' > /dev/null 2>&1 &

### functional enrichment analysis of network results
nohup bash -c '09_functional-enrichment-analysis/03_network-analysis/network-functional-enrichment.R' > /dev/null 2>&1 &

### combine network enrichment plots
nohup bash -c 'Rscript 09_functional-enrichment-analysis/03_network-analysis/03_figures/combine-network-enrichment-plots.R' > /dev/null 2>&1 &