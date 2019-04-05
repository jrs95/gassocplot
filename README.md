# gassocplot
Regional association plots for genetic and epigenetic data.

## Functions
* assoc_plot - plots a regional association plot or fine-mapping probability plot for a single trait within a genomic region based on GRCh37 (hg19) cooridnates.  
* assoc_plot_save - saves a PNG of the assoc_plot.  
* stack_assoc_plot - plots a stacked regional association plot for multiple traits within a genomic region based on GRCh37 (hg19) cooridnates.  
* stack_assoc_plot_save - saves a PNG of the stack_assoc_plot.  

## Installation
1. install.packages("devtools")
2. library(devtools) 
3. install_github("jrs95/gassocplot")
4. library(gassocplot)

## Examples
\#\#\# assoc_plot  
markers <- gassocplot::test_assoc_plot  
head(markers)  
corr <- gassocplot::test_corr # this is correlation not correlation squared  
plot <- assoc_plot(markers, corr)   
assoc_plot_save(plot, "assoc_plot_test.png")  

\#\#\# stack_assoc_plot  
markers <- gassocplot::test_stack_assoc_plot_markers  
head(markers)  
z <- gassocplot::test_stack_assoc_plot_associations  
head(z)  
corr <- gassocplot::test_corr # this is correlation not correlation squared  
plot <- stack_assoc_plot(markers, z, corr, traits=c("Trait 1", "Trait 2"))  
stack_assoc_plot_save(plot, "stack_assoc_plot_test.png", 2)

## Citation
Please cite this R package using the link: https://github.com/jrs95/gassocplot

## Plots

### Regional association plot
![](https://raw.githubusercontent.com/jrs95/utilities/master/assoc_plot_test.png?raw=true)

### Stacked regional association plot
![](https://raw.githubusercontent.com/jrs95/utilities/master/stack_assoc_plot_test.png?raw=true)
