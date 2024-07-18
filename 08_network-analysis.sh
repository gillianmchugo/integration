# network analysis
## network analysis of de results
### genecards
#### search https://www.genecards.org for genes relating to trypano*
#### select export as excel from the export dropdown menu
#### filter to select genes with a relevance score greater than or equal to 1.75
#### convert gene symbols to ensembl ids
Rscript 08_network-analysis/genecards-ids.R

##### found 1036 results, filtered 417 had a relevance score greater than or equal to 1.75 and a valid bovine ensembl id

### innatedb
#### select network analysis from the data analysis dropdown menu on https://www.innatedb.com
#### click upload a file, select prepared genecards results, click open and continue
#### check the box to include interactions predicted by orthology
#### click column 1 and select this column is cross-reference id from the dropdown menu
#### select ensembl as the cross-reference database from the dropdown menu and click ok
#### click next
#### select sif and xggml from the download dropdown menu
#### prepare base network ids
Rscript 08_network-analysis/base-network-ids.R
##### 14077 interactions

### cytoscape
#### open cytoscape and click file > import > network from file... and select sif file from innatedb and click open
#### in the network tab right click on the network collection and click rename network collection..., enter new name and click ok
#### repeat for the network
#### click tools > analyze network and ok
#### in the style tab change the border paint colour to 64AE96 in the rgb tab
#### set the border width to 7.5
#### set the fill colour to #64AE96
#### map the label size to degree with the minimum set as 12 and maximum set as 120
#### set the shape to ellipse and click apply
#### tick the box to lock node width and height
#### map the size to degree with continuous mapping
#### set the minimum node size to 50 and maximum to the maximum degree
#### set the stroke colour to #666666 and transparency to 50
#### click layout > yfiles remove overlaps
#### click file > import > table from file... and select de gene results and click open
#### select to a network collection from the dropdown menu and click on symbol column and the key button and click ok
#### click file > new network > clone current network
#### rename the network collection with the contrast
#### in the style tab click the menu button and copy style..., enter the name of the contrast and click ok
#### map the colour to the coef for the contrast with continuous mapping and set the negative values to E35E2C and positive values to 8E509A
#### map the border colour to the results for the contrast with discrete mapping and set -1 and 1 as the colour for the tissue
#### blood results to D62221
#### liver results to 80B973
#### lymph node results to 4F74C1
#### spleen results to D4B13F
#### in the jactive modules tab select the target network from the dropdown menu
#### click on the adjusted p-value for the contrast to highlight it and click search
#### for each of the resulting modules click tools > analyze network and ok

### identify shared genes between modules
Rscript 08_network-analysis/02_output/shared-module-genes.R

### crop base network plot
Rscript 08_network-analysis/03_figures/01_base-network/crop-base-network.R

### add legends to modules
Rscript 08_network-analysis/03_figures/02_modules/module-legends.R

### combine modules
Rscript 08_network-analysis/03_figures/02_modules/combine-modules.R

### generate bar and upset plots for modules
Rscript 08_network-analysis/03_figures/module-bar-upset.R