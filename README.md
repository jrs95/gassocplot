# gassocplot
Regional association plots for genetic and epigenetic data.

# Functions
* assoc_plot - plots a regional association plot or fine-mapping probability plot for a single trait within a genomic region.  
* assoc_plot_save - saves a PNG of the assoc_plot.  
* stack_assoc_plot - plots a stacked regional association plot for multiple traits within a genomic region.  
* stack_assoc_plot_save - saves a PNG of the stack_assoc_plot.  

# Installation
1. install.packages("devtools")
2. library(devtools) 
3. install_github("jrs95/gassocplot")
4. library(gassocplot)

# Examples
\#\#\# assoc_plot  
markers <- gassocplot::test_assoc_plot  
head(markers)  
corr <- gassocplot::test_corr   
assoc_plot(markers, corr)   

\#\#\# stack_assoc_plot  
markers <- gassocplot::test_stack_assoc_plot_markers  
head(markers)  
z <- gassocplot::test_stack_assoc_plot_associations  
head(z)  
corr <- gassocplot::test_corr   
stack_assoc_plot(markers, z, corr, c("T1", "T2"))   
