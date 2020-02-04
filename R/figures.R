########################################################################
## Regional association plots                                         ##
##                                                                    ##
## James Staley                                                       ##
## University of Bristol                                              ## 
## Email: james.staley@bristol.ac.uk                                  ##
##                                                                    ##
## 19/02/19                                                           ##
########################################################################

##########################################################
##### Recombination plot #####
##########################################################

#' plot_recombination_rate
#'
#' plot_recombination_rate plots the recombination rate against location
#' @param chr chromosome
#' @param x.min start of region
#' @param x.max end of region
#' @import ggplot2
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_recombination_rate <- function(chr, x.min, x.max) {
  recombination.data <- genetic_map[(genetic_map$chr==chr & genetic_map$pos>=x.min & genetic_map$pos<=x.max), ]
  recomb.df <- data.frame(coordinates=recombination.data$pos, y=recombination.data$combined_rate, panel="Recombination Rate", stringsAsFactors=F)
  cols <- c("Recomb. rate"="black","Gene"="#FF3D14","Exon"="#66A300")
  recomb.plot <- ggplot(data = recomb.df, aes(x=coordinates, y=y, colour="black")) + 
  geom_line(aes(colour="black"), colour="steelblue1") + scale_y_continuous(breaks=c(0,25,50,75,100), limits=c(0,100)) + xlab(NULL) + scale_x_continuous(limits=c(x.min,x.max), breaks=NULL) + ylab("Recomb. Rate") + theme_bw() + theme(axis.title.y=element_text(vjust=1.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.title=element_text(size=16), axis.text=element_text(size=14)) + theme(panel.background = element_rect(fill = NA))
  return(recomb.plot)
}

#' plot_recombination_rate_stack
#'
#' plot_recombination_rate_stack plots the recombination rate against location for stacked plots
#' @param chr chromosome
#' @param x.min start of region
#' @param x.max end of region
#' @import ggplot2
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_recombination_rate_stack <- function(chr, x.min, x.max) {
  recombination.data <- genetic_map[(genetic_map$chr==chr & genetic_map$pos>=x.min & genetic_map$pos<=x.max), ]
  recomb.df <- data.frame(coordinates=recombination.data$pos, y=recombination.data$combined_rate, panel="Recombination Rate", stringsAsFactors=F)
  cols <- c("Recomb. rate"="black","Gene"="#FF3D14","Exon"="#66A300")
  recomb.plot <- ggplot(data = recomb.df, aes(x=coordinates, y=y, colour="black")) + 
  geom_line(aes(colour="black"), colour="steelblue1") + scale_y_continuous(breaks=c(0,25,50,75,100), limits=c(0,100)) + xlab(NULL) + scale_x_continuous(limits=c(x.min,x.max), breaks=NULL) + ylab("Recomb. Rate") + theme_bw() + theme(axis.title.y=element_text(vjust=1.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.title=element_text(size=14), axis.text=element_text(size=12)) + theme(panel.background = element_rect(fill = NA))
  return(recomb.plot)
}

##########################################################
##### Gene bar #####
##########################################################

#' plot_gene_zero
#'
#' plot_gene_zero plots a blank gene bar if there are no genes located in the specified region 
#' @param chr chromosome
#' @param x.min start of region
#' @param x.max end of region
#' @import ggplot2
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_gene_zero <- function(chr, x.min, x.max, stack=FALSE){
  genes.df.pos <- data.frame(pos=c(x.min,x.max), y=c(10,5), stringsAsFactors=F, stack=FALSE) 
  plot.genes <- ggplot(data=genes.df.pos, aes(x=pos, y=y)) + theme_bw() + xlab(paste0("Position on chromosome ", chr)) + ylab(" ") +  scale_y_continuous(limits=c(-5,17), breaks=c(8,16), labels=c("      ", "      ")) + theme(axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-0.5), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(limits=c(x.min, x.max)) 
  if(stack==TRUE){
    plot.genes <- plot.genes + theme(axis.title=element_text(size=14), axis.text=element_text(size=12))
  }else{
    plot.genes <- plot.genes + theme(axis.title=element_text(size=16), axis.text=element_text(size=14))
  }
  return(plot.genes)
}

#' plot_gene_two
#'
#' plot_gene_two plots a gene bar with two rows 
#' @param gene.region the gene.region dataset
#' @param chr chromosome
#' @param x.min start of region
#' @param x.max end of region
#' @import ggplot2
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_gene_two <- function(gene.region, chr, x.min, x.max, stack=FALSE) {
  small.gene <- (as.numeric(gene.region$end) - as.numeric(gene.region$start)) < (x.max-x.min)/190
  gene.region$mid.point <- as.numeric(gene.region$start)+(as.numeric(gene.region$end) - as.numeric(gene.region$start))/2
  gene.region$start[small.gene] <- gene.region$mid.point[small.gene] - (x.max-x.min)/380
  gene.region$end[small.gene] <- gene.region$mid.point[small.gene] + (x.max-x.min)/380
  genes.start.stop <- as.matrix(gene.region[, c("start", "end")])  
  genes.df.pos <- data.frame(name=paste0("gene",rep(1:nrow(genes.start.stop), each=2)), pos=as.vector(t(genes.start.stop)), y=(16 - 8*rep(rep(1:2, each=2), ceiling(nrow(genes.start.stop)/2))[1:(2*nrow(genes.start.stop))]), stringsAsFactors=F)
  # plot.pos <- ggplot(data=genes.df.pos, aes(x=pos, y=y)) + geom_point(data=genes.df.pos, mapping=aes(x=pos, y=y, group=name), color="blue4", size=0.5, shape=15) + theme_bw() + xlab(paste0("Position on chromosome ", chr)) + ylab(" ") +  scale_y_continuous(limits=c(-5,17), breaks=c(8,16), labels=c("      ", "      ")) + theme(axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-0.5), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(limits=c(x.min, x.max))
  plot.pos <- ggplot(data=genes.df.pos, aes(x=pos, y=y)) + theme_bw() + xlab(paste0("Position on chromosome ", chr)) + ylab(" ") +  scale_y_continuous(limits=c(-5,17), breaks=c(8,16), labels=c("      ", "      ")) + theme(axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-0.5), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(limits=c(x.min, x.max))
  if(stack==TRUE){
    plot.pos <- plot.pos + theme(axis.title=element_text(size=14), axis.text=element_text(size=12))
  }else{
    plot.pos <- plot.pos + theme(axis.title=element_text(size=16), axis.text=element_text(size=14))
  }
  plot.genes <- plot.pos + geom_line(data=genes.df.pos, aes(x=pos, y=y, group=name), colour="blue4", size=0.8) 
  genes.df.mid.point <- data.frame(name=gene.region$gene, x=as.numeric(gene.region$mid.point), y=(16 - 8*rep(rep(1:2, each=1), ceiling(nrow(gene.region)/2))[1:nrow(gene.region)] + 3.6), stringsAsFactors=F)
  plot.genes <- plot.genes + geom_text(data=genes.df.mid.point, mapping=aes(x=x, y=y, label=name), color="black", size=4, fontface=3) 
  return(plot.genes)
}

#' plot_gene_five
#'
#' plot_gene_five plots a gene bar with five rows 
#' @param gene.region the gene.region dataset
#' @param chr chromosome
#' @param x.min start of region
#' @param x.max end of region
#' @import ggplot2
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_gene_five <- function(gene.region, chr, x.min, x.max, stack=FALSE) {
  small.gene <- (as.numeric(gene.region$end) - as.numeric(gene.region$start)) < (x.max-x.min)/190
  gene.region$mid.point <- as.numeric(gene.region$start)+(as.numeric(gene.region$end) - as.numeric(gene.region$start))/2
  gene.region$start[small.gene] <- gene.region$mid.point[small.gene] - (x.max-x.min)/380
  gene.region$end[small.gene] <- gene.region$mid.point[small.gene] + (x.max-x.min)/380
  genes.start.stop <- as.matrix(gene.region[, c("start", "end")])
  genes.df.pos <- data.frame(name=paste0("gene",rep(1:nrow(genes.start.stop), each=2)), pos=as.vector(t(genes.start.stop)), y=(40 - 8*rep(rep(1:5, each=2), ceiling(nrow(genes.start.stop)/5))[1:(2*nrow(genes.start.stop))]), stringsAsFactors=F)
  # plot.pos <- ggplot(data=genes.df.pos, aes(x=pos, y=y)) + geom_point(data=genes.df.pos, mapping=aes(x=pos, y=y, group=name), color="blue4", size=0.5, shape=15) + theme_bw() + xlab(paste0("Position on chromosome ", chr)) + ylab(" ") +  scale_y_continuous(limits=c(-1,41), breaks=c(8,16), labels=c("      ", "      ")) + theme(axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-0.5), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(limits=c(x.min, x.max))
  plot.pos <- ggplot(data=genes.df.pos, aes(x=pos, y=y)) + theme_bw() + xlab(paste0("Position on chromosome ", chr)) + ylab(" ") +  scale_y_continuous(limits=c(-1,41), breaks=c(8,16), labels=c("      ", "      ")) + theme(axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-0.5), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(limits=c(x.min, x.max))
  if(stack==TRUE){
    plot.pos <- plot.pos + theme(axis.title=element_text(size=14), axis.text=element_text(size=12))
  }else{
    plot.pos <- plot.pos + theme(axis.title=element_text(size=16), axis.text=element_text(size=14))
  }
  plot.genes <- plot.pos + geom_line(data=genes.df.pos, aes(x=pos, y=y, group=name), colour="blue4", size=0.8) 
  genes.df.mid.point <- data.frame(name=gene.region$gene, x=as.numeric(gene.region$mid.point), y=(40 - 8*rep(rep(1:5, each=1), ceiling(nrow(gene.region)/5))[1:nrow(gene.region)] + 3.7), stringsAsFactors=F)
  plot.genes <- plot.genes + geom_text(data=genes.df.mid.point, mapping=aes(x=x, y=y, label=name), color="black", size=3.5, fontface=3) 
  return(plot.genes)
}

#' plot_gene_ten
#'
#' plot_gene_ten plots a gene bar with ten rows 
#' @param gene.region the gene.region dataset
#' @param chr chromosome
#' @param x.min start of region
#' @param x.max end of region
#' @import ggplot2
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_gene_ten <- function(gene.region, chr, x.min, x.max, stack=FALSE) {
  small.gene <- (as.numeric(gene.region$end) - as.numeric(gene.region$start)) < (x.max-x.min)/190
  gene.region$mid.point <- as.numeric(gene.region$start)+(as.numeric(gene.region$end) - as.numeric(gene.region$start))/2
  gene.region$start[small.gene] <- gene.region$mid.point[small.gene] - (x.max-x.min)/380
  gene.region$end[small.gene] <- gene.region$mid.point[small.gene] + (x.max-x.min)/380
  genes.start.stop <- as.matrix(gene.region[, c("start", "end")])
  genes.df.pos <- data.frame(name=paste0("gene",rep(1:nrow(genes.start.stop), each=2)), pos=as.vector(t(genes.start.stop)), y=(80 - 8*rep(rep(1:10, each=2), ceiling(nrow(genes.start.stop)/10))[1:(2*nrow(genes.start.stop))]), stringsAsFactors=F)
  # plot.pos <- ggplot(data=genes.df.pos, aes(x=pos, y=y)) + geom_point(data=genes.df.pos, mapping=aes(x=pos, y=y, group=name), color="blue4", size=0.5, shape=15) + theme_bw() + xlab(paste0("Position on chromosome ", chr)) + ylab(" ") +  scale_y_continuous(limits=c(-1,81), breaks=c(10,20), labels=c("      ", "      ")) + theme(axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-0.5), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(limits=c(x.min, x.max))
  plot.pos <- ggplot(data=genes.df.pos, aes(x=pos, y=y)) + theme_bw() + xlab(paste0("Position on chromosome ", chr)) + ylab(" ") +  scale_y_continuous(limits=c(-1,81), breaks=c(10,20), labels=c("      ", "      ")) + theme(axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-0.5), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(limits=c(x.min, x.max))
  if(stack==TRUE){
    plot.pos <- plot.pos + theme(axis.title=element_text(size=14), axis.text=element_text(size=12))
  }else{
    plot.pos <- plot.pos + theme(axis.title=element_text(size=16), axis.text=element_text(size=14))
  }
  plot.genes <- plot.pos + geom_line(data=genes.df.pos, aes(x=pos, y=y, group=name), colour="blue4", size=0.8) 
  genes.df.mid.point <- data.frame(name=gene.region$gene, x=as.numeric(gene.region$mid.point), y=(80 - 8*rep(rep(1:10, each=1), ceiling(nrow(gene.region)/10))[1:nrow(gene.region)] + 3.5), stringsAsFactors=F)
  plot.genes <- plot.genes + geom_text(data=genes.df.mid.point, mapping=aes(x=x, y=y, label=name), color="black", size=3, fontface=3) 
  return(plot.genes)
}

#' plot_gene_fifteen
#'
#' plot_gene_fifteen plots a gene bar with fifteen rows 
#' @param gene.region the gene.region dataset
#' @param chr chromosome
#' @param x.min start of region
#' @param x.max end of region
#' @import ggplot2
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_gene_fifteen <- function(gene.region, chr, x.min, x.max, stack=FALSE) {
  small.gene <- (as.numeric(gene.region$end) - as.numeric(gene.region$start)) < (x.max-x.min)/190
  gene.region$mid.point <- as.numeric(gene.region$start)+(as.numeric(gene.region$end) - as.numeric(gene.region$start))/2
  gene.region$start[small.gene] <- gene.region$mid.point[small.gene] - (x.max-x.min)/380
  gene.region$end[small.gene] <- gene.region$mid.point[small.gene] + (x.max-x.min)/380
  genes.start.stop <- as.matrix(gene.region[, c("start", "end")])
  genes.df.pos <- data.frame(name=paste0("gene",rep(1:nrow(genes.start.stop), each=2)), pos=as.vector(t(genes.start.stop)), y=(120 - 8*rep(rep(1:15, each=2), ceiling(nrow(genes.start.stop)/15))[1:(2*nrow(genes.start.stop))]), stringsAsFactors=F)
  # plot.pos <- ggplot(data=genes.df.pos, aes(x=pos, y=y)) + geom_point(data=genes.df.pos, mapping=aes(x=pos, y=y, group=name), color="blue4", size=0.3, shape=15) + theme_bw() + xlab(paste0("Position on chromosome ", chr)) + ylab(" ") +  scale_y_continuous(limits=c(-1,121), breaks=c(10,20), labels=c("      ", "      ")) + theme(axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-0.5), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(limits=c(x.min, x.max))
  plot.pos <- ggplot(data=genes.df.pos, aes(x=pos, y=y)) + theme_bw() + xlab(paste0("Position on chromosome ", chr)) + ylab(" ") +  scale_y_continuous(limits=c(-1,121), breaks=c(10,20), labels=c("      ", "      ")) + theme(axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-0.5), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(limits=c(x.min, x.max))
  if(stack==TRUE){
    plot.pos <- plot.pos + theme(axis.title=element_text(size=14), axis.text=element_text(size=12))
  }else{
    plot.pos <- plot.pos + theme(axis.title=element_text(size=16), axis.text=element_text(size=14))
  }
  plot.genes <- plot.pos + geom_line(data=genes.df.pos, aes(x=pos, y=y, group=name), colour="blue4", size=0.7) 
  genes.df.mid.point <- data.frame(name=gene.region$gene, x=as.numeric(gene.region$mid.point), y=(120 - 8*rep(rep(1:15, each=1), ceiling(nrow(gene.region)/10))[1:nrow(gene.region)] + 3.5), stringsAsFactors=F)
  plot.genes <- plot.genes + geom_text(data=genes.df.mid.point, mapping=aes(x=x, y=y, label=name), color="black", size=2, fontface=3) 
  return(plot.genes)
}

##########################################################
##### Assoc plot #####
##########################################################

#' plot_assoc
#'
#' plot_assoc plots a scatter graph of associations (e.g. log10 p-values) 
#' @param data data.frame with markername (marker), chromosome (chr), position (pos) and association statistics (stats)
#' @param corr correlation matrix between markers
#' @param corr.top correlation statistics between the top marker and the rest of the markers
#' @param x.min start of region
#' @param x.max end of region
#' @param top.marker the top associated marker, i.e. the marker with the largest -log10p or probability
#' @param ylab the y-axis label
#' @param type the type of the plot either log10p or probabilities
#' @import ggplot2
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_assoc <- function(data, corr=NULL, corr.top=NULL, x.min, x.max, top.marker=NULL, ylab, type="log10p"){
  if(is.null(corr) & is.null(corr.top)) stop("no correlation statistics were input")
  if(is.null(corr) & !is.null(corr.top) & is.null(top.marker)) stop("top.marker must be defined if corr.top is provided")
  miss <- is.na(data$stats)
  if(!is.null(corr)){corr <- corr[!miss, !miss]}
  if(!is.null(corr.top)){corr.top <- corr.top[!miss]}  
  data <- data[!miss,]
  if(length(top.marker)!=0){
    top_marker <- which(data$marker==top.marker)
    if(length(top_marker)>1){top_marker <- sample(top_marker, 1); if(is.null(corr) & !is.null(corr.top)) warning("top.marker maps to multiple markers")}
    if(length(top_marker)==0){top_marker <- max.col(t(data$stats))} 
    lead_marker <- data[top_marker,]  
    ov_lead_marker <- data[max.col(t(data$stats)),]
    if((lead_marker$stats/ov_lead_marker$stats)>0.975){geomtext <- T}else{geomtext <- F}
    lead_marker$label_pos <- lead_marker$pos
    if((x.max-lead_marker$pos)<10000){lead_marker$label_pos <- lead_marker$pos - 0.025*(x.max-x.min)}
    if((lead_marker$pos-x.min)<10000){lead_marker$label_pos <- lead_marker$pos + 0.025*(x.max-x.min)}
  }else{ 
    top_marker <- max.col(t(data$stats))
    lead_marker <- data[top_marker,]  
    lead_marker$label_pos <- lead_marker$pos
    if((x.max-lead_marker$pos)<10000){lead_marker$label_pos <- lead_marker$pos - 0.025*(x.max-x.min)}
    if((lead_marker$pos-x.min)<10000){lead_marker$label_pos <- lead_marker$pos + 0.025*(x.max-x.min)}
    geomtext <- T
  }
  if(!is.null(corr)){r2 <- corr[,top_marker]^2}else{r2 <- corr.top^2}
  data$r2 <- "miss"
  data$r2[r2>=0 & r2<0.2 & !is.na(r2)] <- "0.0-0.2"
  data$r2[r2>=0.2 & r2<0.4 & !is.na(r2)] <- "0.2-0.4"
  data$r2[r2>=0.4 & r2<0.6 & !is.na(r2)] <- "0.4-0.6"
  data$r2[r2>=0.6 & r2<0.8 & !is.na(r2)] <- "0.6-0.8"
  data$r2[r2>=0.8 & r2<=1 & !is.na(r2)] <- "0.8-1.0" 
  data$r2 <- factor(data$r2, levels=c("miss", "0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"))
  ylim <- max((max(data$stats)+0.1*max(data$stats)),1)
  marker.plot <- ggplot(aes(x=pos, y=stats), data=data) + geom_point(aes(fill=r2), pch=21, size=3.5) + scale_fill_manual(values=c("#DCDCDC", "#66FFFF", "#66FF66", "#FFCC00", "#FF9933", "#CC3300"), drop=FALSE) + geom_point(data=lead_marker, aes(pos,stats), pch=23, colour="black", fill="purple", size=4)  + theme_bw() + ylab(ylab) + xlab(NULL) + scale_y_continuous(limits=c(0,ylim)) + theme(axis.title.y=element_text(vjust=2.25, size=16), axis.text=element_text(size=14)) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + scale_x_continuous(limits=c(x.min,x.max), breaks=NULL) + theme(axis.title=element_text(size=10)) + theme(legend.text=element_text(size=11), legend.title=element_text(size=12), legend.background = element_rect(colour = "black")) + theme(panel.background=element_rect(fill=NA)) + theme(legend.position="bottom") + guides(fill = guide_legend(nrow = 1))
  if(geomtext){marker.plot <- marker.plot + geom_text(data=lead_marker, aes(x=label_pos,y=stats,label=marker), vjust=-1, hjust=0.5, size=4.5)}else{if(lead_marker$stats[1]/ylim>=0.3){marker.plot <- marker.plot + geom_label(data=lead_marker, aes(x=label_pos,y=stats,label=marker), label.r=unit(0, "lines"), nudge_y=(-0.05*ylim), size=4.5, alpha=1)}else{marker.plot <- marker.plot + geom_label(data=lead_marker, aes(x=label_pos,y=stats,label=marker), label.r=unit(0, "lines"), nudge_y=(0.05*ylim), size=4.5, alpha=1)}}
  if(type=="prob"){suppressMessages(marker.plot <- marker.plot + scale_y_continuous(limits=c(0,ylim), breaks=c(0, 0.25, 0.5, 0.75, 1)))}
  return(marker.plot)
}

##########################################################
##### Legend #####
##########################################################

#' g_legend
#'
#' g_legend 
#' @param gplot a ggplot
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
g_legend<-function(gplot){
  tmp <- ggplot_gtable(ggplot_build(gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

##########################################################
##### Combined plot #####
##########################################################

#' plot_assoc_combined
#'
#' plot_assoc_combined combines the recombination plot, gene bar and association scatter graph  
#' @param recombination.plot recombination plot 
#' @param gene.plot gene bar
#' @param marker.plot association scatter plot
#' @param title title of the plot
#' @param subtitle subtitle of the plot
#' @param ngenes number of genes in the genomic region
#' @param r2_legend add r2 legend
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_assoc_combined <- function(recombination.plot, gene.plot, marker.plot, title=NULL, subtitle=NULL, ngenes, r2_legend=TRUE){
  legend <- g_legend(marker.plot); marker.plot <- marker.plot + theme(legend.position="none")
  g1 <- ggplot_gtable(ggplot_build(recombination.plot))
  g2 <- ggplot_gtable(ggplot_build(marker.plot))
  pp <- c(subset(g1$layout, name == "panel", se = t:r))
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)
  ia <- which(g1$layout$name == "axis-l")
  ga <- g1$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)
  ax$grobs <- rev(ax$grobs)
  ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
  g <- gtable_add_cols(g, g1$widths[g1$layout[ia, ]$l], length(g$widths) - 1)
  g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
  g$grobs[[3]] <- g2$grobs[[3]]
  g$grobs[[13]] <- g2$grobs[[13]]
  g <- gtable_add_cols(g, g1$widths[g1$layout[ia, ]$l], length(g$widths) - 1)
  if(ngenes<=10){
    g <- gtable_add_grob(g, list(textGrob("Recombination Rate (cM/Mb)", rot = -90, gp = gpar(col="black", fontsize=16))), pp$t, length(g$widths) - 1, pp$b)
  }else{
    g <- gtable_add_grob(g, list(textGrob("Recombination Rate", rot = -90, gp = gpar(col="black", fontsize=16))), pp$t, length(g$widths) - 1, pp$b)
  }
  g3 <- ggplot_gtable(ggplot_build(gene.plot))
  g3 <- gtable_add_grob(g3, g3$grobs[[which(g3$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)
  ia <- which(g3$layout$name == "axis-l")
  ga <- g3$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)
  ax$grobs <- rev(ax$grobs)
  g3 <- gtable_add_cols(g3, g3$widths[g3$layout[ia, ]$l], length(g3$widths) - 1)
  g3 <- gtable_add_grob(g3, ax, pp$t, length(g3$widths) - 1, pp$b)
  g3 <- gtable_add_cols(g3, g3$widths[g3$layout[ia, ]$l], length(g3$widths) - 1)
  g3 <- gtable_add_grob(g3, list(textGrob("", rot = -90, gp = gpar(fontsize=16, col = gray(.88)))), pp$t, length(g3$widths) - 1, pp$b)
  g <- gtable:::rbind_gtable(g, g3, "last")
  panels <- g$layout$t[grep("panel", g$layout$name)]
  g$heights[panels[1]] <- unit(3,"null") 
  if(ngenes<=5){g$heights[panels[2]] <- unit(0.5,"null")}
  if(ngenes>5 & ngenes<=10){g$heights[panels[2]] <- unit(1.2,"null")}
  if(ngenes>10){g$heights[panels[2]] <- unit(2,"null")}
  if(!is.null(subtitle)){
    gt1 <- textGrob(subtitle,gp=gpar(fontsize=16))
    g <- gtable_add_rows(g, heights = grobHeight(gt1)*2.5, pos = 0)
    g <- gtable_add_grob(g, gt1, 1, 1, 1, ncol(g))
  }
  if(!is.null(title)){
    gt <- textGrob(title,gp=gpar(fontsize=20, fontface="italic"))
    g <- gtable_add_rows(g, heights = grobHeight(gt)*1.5, pos = 0)
    g <- gtable_add_grob(g, gt, 1, 1, 1, ncol(g))
  }
  g <- gtable_add_padding(g, unit(0.3, "cm"))
  if(r2_legend==TRUE){
    lheight <- sum(legend$height)*1.5
    g <- grid.arrange(g, legend, ncol = 1, heights = unit.c(unit(1, "npc") - lheight, lheight))
  }
  return(g)
}

##########################################################
##### Assoc plot #####
##########################################################

#' assoc_plot
#'
#' assoc_plot plots a scatter graph of associations (e.g. log10 p-values)
#' @param data data.frame with markername (marker), chromosome (chr), position (pos) and either z-statistics (z) or probabilities (prob) columns
#' @param corr correlation matrix between markers
#' @param corr.top correlation statistics between the top marker and the rest of the markers
#' @param ylab the y-axis label
#' @param title title of the plot
#' @param subtitle subtitle of the plot
#' @param type the type of the plot either log10p or probabilities
#' @param x.min start of region
#' @param x.max end of region
#' @param top.marker the top associated marker, i.e. the marker with the largest -log10p or probability
#' @param legend add r2 legend
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
assoc_plot <- function(data, corr=NULL, corr.top=NULL, ylab=NULL, title=NULL, subtitle=NULL, type="log10p", x.min=NULL, x.max=NULL, top.marker=NULL, legend=TRUE){
  
  # Error messages
  if(!(type=="log10p" | type=="prob")) stop("the type of plot has to be either log10p or prob")
  if(type=="log10p"){
    if(any(names(data)!=c("marker", "chr", "pos", "z"))) stop("dataset needs to include marker, chr, pos and z columns in that order")
  }else{
    if(any(names(data)!=c("marker", "chr", "pos", "prob"))) stop("dataset needs to include marker, chr, pos and prob columns in that order")
  }
  if(!is.null(corr)){if(ncol(corr)!=nrow(data) | nrow(corr)!=nrow(data)) stop("corr has to have the same dimensions as the number of rows in the markers dataset")}
  # if(any(rownames(corr)!=data$marker)) stop("corr has to have the same markers in the same order as the dataset")
  if(length(unique(data$chr))>1) stop("there should only be markers from one chromosome in the markers dataset") 
  if(!(data$chr[1] %in% 1:22)) stop("the plotting tool is only for autosomal chromosomes") 
  if(any(is.na(data))) stop("there are missing values in the dataset") 
  if(class(data$pos)!="integer") stop("the pos variable has to be an integer")
  if(is.null(corr) & !is.null(corr.top) & is.null(top.marker)) stop("top.marker must be defined if corr.top is provided")
  if(is.null(corr) & !is.null(corr.top)){if(length(corr.top)!=nrow(data)) stop("corr.top has to have the same length as the number of rows in the markers dataset")}
  if(!is.null(top.marker) & length(which(top.marker==data$marker))==0) stop("top.marker is not contained in the markers dataset")
  if(!is.null(top.marker) & length(which(top.marker==data$marker))>1) stop("top.marker maps to multiple markers in the markers dataset")

  # Dataset
  if(type=="log10p"){
    mlog10p <- -(log(2) + pnorm(-abs(data$z), log.p=T))/log(10)
    mlog10p[mlog10p>1000] <- 1000
    data$stats <- mlog10p
  }else{data$stats <- data$prob}
  data <- data[,c("marker", "chr", "pos", "stats")]
  data$marker <- as.character(data$marker)
  chr <- as.integer(data$chr[1])
  if(is.null(x.min)){x.min <- min(as.integer(data$pos))}
  if(is.null(x.max)){x.max <- max(as.integer(data$pos))}
  if((x.max - x.min)>10000000) stop("the plotting tool can plot a maximum of 10MB")

  # Genes
  gene.region <- gassocplot::genes[gassocplot::genes$chr==chr & !(gassocplot::genes$end<x.min) & !(gassocplot::genes$start>x.max),]
  gene.region$start[gene.region$start<x.min] <- x.min
  gene.region$end[gene.region$end>x.max] <- x.max
  gene.region <- gene.region[with(gene.region, order(start)), ]
  ngenes <- nrow(gene.region)

  # Max and min
  x.min <- x.min - 0.02*(x.max - x.min)
  x.max <- x.max + 0.02*(x.max - x.min)
  
  # Correlation matrix
  if(is.null(corr) & is.null(corr.top)){r2_legend <- FALSE; corr <- matrix(NA, nrow=nrow(markers), ncol=nrow(markers))}
  
  # Recombination plot
  recombination.plot <- plot_recombination_rate(chr, x.min, x.max)

  # Gene plot
  if(ngenes==0){gene.plot <- plot_gene_zero(chr, x.min, x.max)}
  if(ngenes>0 & ngenes<=5){gene.plot <- plot_gene_two(gene.region, chr, x.min, x.max)}
  if(ngenes>5 & ngenes<=10){gene.plot <- plot_gene_five(gene.region, chr, x.min, x.max)}
  if(ngenes>10 & ngenes<=25){gene.plot <- plot_gene_ten(gene.region, chr, x.min, x.max)}
  if(ngenes>25){gene.plot <- plot_gene_fifteen(gene.region, chr, x.min, x.max)}
  
  # Marker plot
  data$chr <- as.integer(data$chr)
  data$pos <- as.integer(data$pos)
  if(type=="log10p"){ylab <- expression("-log"["10"]*paste("(",italic("p"),")"))}else{if(is.null(ylab)){ylab <- "Probability"}}  
  marker.plot <- plot_assoc(data, corr, corr.top, x.min, x.max, top.marker, ylab, type)
  
  # Combined plot
  combined.plot <- plot_assoc_combined(recombination.plot, gene.plot, marker.plot, title, subtitle, ngenes, legend)

  return(combined.plot)
}

##########################################################
##### Save assoc plot #####
##########################################################

#' assoc_plot_save
#'
#' assoc_plot_save saves a png of the assoc_plot with the correct dimensions
#' @param x the plot
#' @param file the filepath
#' @param width the width of the plot
#' @param height the height of the plot
#' @param dpi the resolution of the plot
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
assoc_plot_save <- function(x, file, width=9, height=7, dpi=500){
  suppressGraphics(ggsave(file, plot=grid.draw(x), width=width, height=height, units="in", limitsize=F, dpi=dpi))
}

##########################################################
##### Stack assoc plot #####
##########################################################

#' plot_assoc_stack
#'
#' plot_assoc_stack plots a scatter graph of -log10p against chromosmal location without gene bar
#' @param data data.frame with markername (marker), chromosome (chr), position (pos) and association statistics (stats)
#' @param corr correlation matrix between markers
#' @param corr.top correlation statistics between the top marker and the rest of the markers
#' @param x.min start of region
#' @param x.max end of region
#' @param top.marker the top associated marker, i.e. the marker with the largest -log10p or probability
#' @param ylab the y-axis label
#' @param type the type of the plot either log10p or probabilities
#' @import ggplot2
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_assoc_stack <- function(data, corr=NULL, corr.top=NULL, x.min, x.max, top.marker=NULL, ylab, type="log10p"){
  if(is.null(corr) & is.null(corr.top)) stop("no correlation statistics were input")
  miss <- is.na(data$stats)
  if(!is.null(corr)){corr <- corr[!miss, !miss]}
  if(!is.null(corr.top)){corr.top <- corr.top[!miss]}  
  data <- data[!miss,]
  if(length(top.marker)!=0){
    top_marker <- which(data$marker==top.marker)
    if(length(top_marker)>1){top_marker <- sample(top_marker, 1); if(is.null(corr) & !is.null(corr.top)) warning("top.marker maps to multiple markers")}
    if(length(top_marker)==0){top_marker <- max.col(t(data$stats))} 
    lead_marker <- data[top_marker,]  
    ov_lead_marker <- data[max.col(t(data$stats)),]
    if((lead_marker$stats/ov_lead_marker$stats)>0.975){geomtext <- T}else{geomtext <- F}
    lead_marker$label_pos <- lead_marker$pos
    if((x.max-lead_marker$pos)<10000){lead_marker$label_pos <- lead_marker$pos - 0.025*(x.max-x.min)}
    if((lead_marker$pos-x.min)<10000){lead_marker$label_pos <- lead_marker$pos + 0.025*(x.max-x.min)}
  }else{ 
    top_marker <- max.col(t(data$stats))
    lead_marker <- data[top_marker,]  
    lead_marker$label_pos <- lead_marker$pos
    if((x.max-lead_marker$pos)<10000){lead_marker$label_pos <- lead_marker$pos - 0.025*(x.max-x.min)}
    if((lead_marker$pos-x.min)<10000){lead_marker$label_pos <- lead_marker$pos + 0.025*(x.max-x.min)}
    geomtext <- T
  }
  if(!is.null(corr)){r2 <- corr[,top_marker]^2}else{r2 <- corr.top^2}
  data$r2 <- "miss"
  data$r2[r2>=0 & r2<0.2 & !is.na(r2)] <- "0.0-0.2"
  data$r2[r2>=0.2 & r2<0.4 & !is.na(r2)] <- "0.2-0.4"
  data$r2[r2>=0.4 & r2<0.6 & !is.na(r2)] <- "0.4-0.6"
  data$r2[r2>=0.6 & r2<0.8 & !is.na(r2)] <- "0.6-0.8"
  data$r2[r2>=0.8 & r2<=1 & !is.na(r2)] <- "0.8-1.0"
  data$r2 <- factor(data$r2, levels=c("miss", "0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"))
  ylim <- max((max(data$stats)+0.2*max(data$stats)),1)
  marker.plot <- ggplot(aes(x=pos,y=stats), data=data) + geom_point(aes(fill=r2), pch=21, size=3) + scale_fill_manual(values=c("#DCDCDC", "#66FFFF", "#66FF66", "#FFCC00", "#FF9933", "#CC3300"), drop=FALSE) + geom_point(data=lead_marker, aes(x=pos,y=stats), pch=23, colour="black", fill="purple", size=4) + theme_bw() +  ylab(ylab) + xlab(NULL) + scale_y_continuous(limits=c(0,ylim)) + theme(axis.title.y=element_text(vjust=2.25, size=14), axis.text=element_text(size=12)) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + scale_x_continuous(limits=c(x.min, x.max), breaks=NULL) + theme(axis.title=element_text(size=10)) + theme(legend.text=element_text(size=10), legend.title=element_text(size=12), legend.background = element_rect(colour = "black")) + theme(panel.background=element_rect(fill=NA)) + theme(legend.position="bottom") + guides(fill = guide_legend(nrow = 1))
  if(geomtext){marker.plot <- marker.plot + geom_text(data=lead_marker, aes(x=label_pos,y=stats,label=marker), vjust=-1, hjust=0.5, size=4.5)}else{if(lead_marker$stats[1]/ylim>=0.3){marker.plot <- marker.plot + geom_label(data=lead_marker, aes(x=label_pos,y=stats,label=marker), label.r=unit(0, "lines"), nudge_y=(-0.1*ylim), size=4.5, alpha=1)}else{marker.plot <- marker.plot + geom_label(data=lead_marker, aes(x=label_pos,y=stats,label=marker), label.r=unit(0, "lines"), nudge_y=(0.1*ylim), size=4.5, alpha=1)}}
  if(type=="prob"){suppressMessages(marker.plot <- marker.plot + scale_y_continuous(limits=c(0,ylim), breaks=c(0, 0.25, 0.5, 0.75, 1)))}
  return(marker.plot)
}

##########################################################
##### Regional association plot (without gene bar) #####
##########################################################

#' plot_regional_assoc
#'
#' plot_regional_assoc plots the regional association plot 
#' @param recombination.plot recombination plot 
#' @param marker.plot association scatter plot
#' @param title title of the plot
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_regional_assoc <- function(recombination.plot, marker.plot, title){
  marker.plot <- marker.plot + theme(legend.position="none")
  g1 <- ggplot_gtable(ggplot_build(recombination.plot))
  g2 <- ggplot_gtable(ggplot_build(marker.plot))
  pp <- c(subset(g1$layout, name == "panel", se = t:r))
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)
  ia <- which(g1$layout$name == "axis-l")
  ga <- g1$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)
  ax$grobs <- rev(ax$grobs)
  ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
  g <- gtable_add_cols(g, g1$widths[g1$layout[ia, ]$l], length(g$widths) - 1)
  g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
  g$grobs[[3]] <- g2$grobs[[3]]
  g$grobs[[13]] <- g2$grobs[[13]]
  g <- gtable_add_cols(g, g1$widths[g1$layout[ia, ]$l], length(g$widths) - 1)
  g <- gtable_add_grob(g, list(textGrob("Recomb. Rate", rot = -90, gp = gpar(col="black", fontsize=14))), pp$t, length(g$widths) - 1, pp$b)
  gt <- textGrob(title,gp=gpar(fontsize=16, fontface="italic"))
  g <- gtable_add_rows(g, heights = grobHeight(gt)*1.5, pos = 0)
  g <- gtable_add_grob(g, gt, 1, 1, 1, ncol(g))
  g <- gtable_add_padding(g, unit(0.3, "cm"))
  return(g)
}

##########################################################
##### Association plot (with gene bar) #####
##########################################################

#' plot_regional_gene_assoc
#'
#' plot_regional_gene_assoc plots the regional association plot 
#' @param recombination.plot recombination plot 
#' @param gene.plot gene bar
#' @param marker.plot association scatter plot
#' @param title title of the plot
#' @param ngenes number of genes in the genomic region
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_regional_gene_assoc <- function(recombination.plot, marker.plot, gene.plot, title, ngenes){
  marker.plot <- marker.plot + theme(legend.position="none")
  g1 <- ggplot_gtable(ggplot_build(recombination.plot))
  g2 <- ggplot_gtable(ggplot_build(marker.plot))
  pp <- c(subset(g1$layout, name == "panel", se = t:r))
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)
  ia <- which(g1$layout$name == "axis-l")
  ga <- g1$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)
  ax$grobs <- rev(ax$grobs)
  ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
  g <- gtable_add_cols(g, g1$widths[g1$layout[ia, ]$l], length(g$widths) - 1)
  g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
  g$grobs[[3]] <- g2$grobs[[3]]
  g$grobs[[13]] <- g2$grobs[[13]]
  g <- gtable_add_cols(g, g1$widths[g1$layout[ia, ]$l], length(g$widths) - 1)
  g <- gtable_add_grob(g, list(textGrob("Recomb. Rate", rot = -90, gp = gpar(col="black", fontsize=14))), pp$t, length(g$widths) - 1, pp$b)
  g3 <- ggplot_gtable(ggplot_build(gene.plot))
  g3 <- gtable_add_grob(g3, g3$grobs[[which(g3$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)
  ia <- which(g3$layout$name == "axis-l")
  ga <- g3$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)
  ax$grobs <- rev(ax$grobs)
  g3 <- gtable_add_cols(g3, g3$widths[g3$layout[ia, ]$l], length(g3$widths) - 1)
  g3 <- gtable_add_grob(g3, ax, pp$t, length(g3$widths) - 1, pp$b)
  g3 <- gtable_add_cols(g3, g3$widths[g3$layout[ia, ]$l], length(g3$widths) - 1)
  g3 <- gtable_add_grob(g3, list(textGrob("", rot = -90, gp = gpar(fontsize=16, col = gray(.88)))), pp$t, length(g3$widths) - 1, pp$b)
  g <- gtable:::rbind_gtable(g, g3, "last")
  panels <- g$layout$t[grep("panel", g$layout$name)]
  g$heights[panels[1]] <- unit(3,"null") 
  if(ngenes<=5){g$heights[panels[2]] <- unit(0.75,"null")}
  if(ngenes>5 & ngenes<=10){g$heights[panels[2]] <- unit(1.8,"null")}
  if(ngenes>10){g$heights[panels[2]] <- unit(2.5,"null")}
  gt <- textGrob(title,gp=gpar(fontsize=16, fontface="italic"))
  g <- gtable_add_rows(g, heights = grobHeight(gt)*1.5, pos = 0)
  g <- gtable_add_grob(g, gt, 1, 1, 1, ncol(g))
  g <- gtable_add_padding(g, unit(0.3, "cm"))
  return(g)
}

##########################################################
##### Add legend to regional association plots #####
##########################################################

#' add_g_legend
#'
#' add_g_legend
#' @param g a ggplot
#' @param legend legend
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
add_g_legend <- function(g, legend){
  lheight <- sum(legend$height)*1.5
  g <- grid.arrange(g, legend, ncol = 1, heights = unit.c(unit(1, "npc") - lheight, lheight))
  return(g)
}

##########################################################
##### Stacked regional association plot #####
##########################################################

#' stack_assoc_plot
#'
#' stack_assoc_plot plots stacked regional association plots
#' @param markers data.frame of markers with markername (marker), chromosome (chr) and position (pos) 
#' @param z matrix of Z-scores or probabilities with one column for each trait
#' @param corr matrix of correlation statistics between markers
#' @param corr.top correlation statistics between the top marker and the rest of the markers
#' @param traits trait names
#' @param ylab the y-axis label
#' @param type the type of the plot either log10p or probabilities
#' @param x.min start of region
#' @param x.max end of region
#' @param top.marker the top associated marker, i.e. the marker with the largest -log10p or probability
#' @param legend add r2 legend
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
stack_assoc_plot <- function(markers, z, corr=NULL, corr.top=NULL, traits, ylab=NULL, type="log10p", x.min=NULL, x.max=NULL, top.marker=NULL, legend=TRUE){
  
  # Error messages
  if(!(type=="log10p" | type=="prob")) stop("the type of plot has to be either log10p or prob")
  if(length(traits)!=ncol(z)) stop("the number of traits is not the same as the number of columns for the Z-scores")
  if(nrow(markers)!=nrow(z)) stop("the number of markers is not the same as the number of rows for the Z-scores")
  if(!is.null(corr)){if(ncol(corr)!=nrow(markers) | nrow(corr)!=nrow(markers)) stop("corr has to have the same dimensions as the number of rows in the markers dataset")}
  if(!is.null(corr.top)){if(length(corr.top)!=nrow(markers)) stop("corr.top has to have the same length as the number of rows in the markers dataset")}
  # if(any(rownames(corr)!=markers$marker)) stop("corr has to have the same markers in the same order as the markers dataset")
  if(any(names(markers)!=c("marker", "chr", "pos"))) stop("dataset needs to include marker, chr and pos columns in that order")
  if(length(unique(markers$chr))>1) stop("there should only be markers from one chromosome in the markers dataset")   
  if(!(markers$chr[1] %in% 1:22)) stop("the plotting tool is only for autosomal chromosomes") 
  if(any(is.na(markers))) stop("there are missing markers in your marker dataset") 
  # if(any(is.na(z))) stop("there are missing values in the Z-score matrix")
  if(class(markers$pos)!="integer") stop("the pos variable has to be an integer")
  if(is.null(corr) & !is.null(corr.top) & is.null(top.marker)) stop("top.marker must be defined if corr.top is provided")
  if(is.null(corr) & !is.null(corr.top)){if(length(corr.top)!=nrow(markers)) stop("corr.top has to have the same length as the number of rows in the markers dataset")}
  if(!is.null(top.marker) & length(which(top.marker==markers$marker))==0) stop("top.marker is not contained in the markers dataset")
  if(!is.null(top.marker) & length(which(top.marker==markers$marker))>1) stop("top.marker maps to multiple markers in the markers dataset")
  
  # Coerce data
  markers$marker <- as.character(markers$marker)
  chr <- as.integer(markers$chr[1])
  r2_legend <- legend
  if(is.null(x.min)){x.min <- min(as.integer(markers$pos))}
  if(is.null(x.max)){x.max <- max(as.integer(markers$pos))}
  if((x.max - x.min)>10000000) stop("the plotting tool can plot a maximum of 10MB")
  
  # mlog10p
  if(type=="log10p"){
    mlog10p <- suppressWarnings(apply(z, 2, function(x){-(log(2) + pnorm(-abs(x), log.p=T))/log(10)}))
    mlog10p[mlog10p>1000 & !is.na(mlog10p)] <- 1000
  }
  if(type=="log10p"){ylab <- expression("-log"["10"]*paste("(",italic("p"),")"))}else{if(is.null(ylab)){ylab <- "Probability"}}
 
  # Genes
  gene.region <- gassocplot::genes[gassocplot::genes$chr==chr & !(gassocplot::genes$end<x.min) & !(gassocplot::genes$start>x.max),]
  gene.region$start[gene.region$start<x.min] <- x.min
  gene.region$end[gene.region$end>x.max] <- x.max
  gene.region <- gene.region[with(gene.region, order(start)), ]
  ngenes <- nrow(gene.region)

  # Max and min
  x.min <- x.min - 0.02*(x.max - x.min)
  x.max <- x.max + 0.02*(x.max - x.min)
  
  # Correlation matrix
  if(is.null(corr) & is.null(corr.top)){r2_legend <- FALSE; corr <- matrix(NA, nrow=nrow(markers), ncol=nrow(markers))}

  # Recombination plot
  recombination.plot <- plot_recombination_rate_stack(chr, x.min, x.max)

  # Gene plot
  if(ngenes==0){gene.plot <- plot_gene_zero(chr, x.min, x.max, stack=TRUE)}
  if(ngenes>0 & ngenes<=5){gene.plot <- plot_gene_two(gene.region, chr, x.min, x.max, stack=TRUE)}
  if(ngenes>5 & ngenes<=10){gene.plot <- plot_gene_five(gene.region, chr, x.min, x.max, stack=TRUE)}
  if(ngenes>10 & ngenes<=25){gene.plot <- plot_gene_ten(gene.region, chr, x.min, x.max, stack=TRUE)}
  if(ngenes>25){gene.plot <- plot_gene_fifteen(gene.region, chr, x.min, x.max, stack=TRUE)}

  # Top marker
  if(length(top.marker)!=0){if(is.na(top.marker)){top.marker <- NULL}}
    
  # Association plot
  for(i in length(traits):1){
    if(type=="log10p"){
      data <- data.frame(marker=markers$marker, chr=as.integer(markers$chr), pos=as.integer(markers$pos), stats=mlog10p[,i], stringsAsFactors=F)
    }else{
      data <- data.frame(marker=markers$marker, chr=as.integer(markers$chr), pos=as.integer(markers$pos), stats=z[,i], stringsAsFactors=F)    
    }
    marker.plot <- plot_assoc_stack(data, corr, corr.top, x.min, x.max, top.marker, ylab, type)
    legend <- g_legend(marker.plot)
    if(i==length(traits)){g <- plot_regional_gene_assoc(recombination.plot, marker.plot, gene.plot, traits[i], ngenes)}
    if(i<length(traits)){
      g1 <- plot_regional_assoc(recombination.plot, marker.plot, traits[i])
      g <- gtable:::rbind_gtable(g1, g, "last")
      panels <- g$layout$t[grep("panel", g$layout$name)]
      g$heights[panels[1]] <- unit(3,"null") 
    } 
  }

  # Combined plot
  if(r2_legend==T){
    combined.plot <- add_g_legend(g, legend)
  }else{
    combined.plot <- g
  }

  return(combined.plot)
}

##########################################################
##### Save assoc plot #####
##########################################################

#' stack_assoc_plot_save
#'
#' stack_assoc_plot_save saves a png of the assoc_plot with the correct dimensions
#' @param x the plot
#' @param file the filepath
#' @param n_traits the filepath
#' @param width the width of the plot
#' @param height the height of the plot
#' @param dpi the resolution of the plot
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
stack_assoc_plot_save <- function(x, file, n_traits, width=NULL, height=NULL, dpi=500){
  if(is.null(width)){width <- 8}
  if(is.null(height)){height <- 3 + 3*n_traits}  
  suppressGraphics(ggsave(file, plot=grid.draw(x), width=width, height=height, units="in", limitsize=F, dpi=dpi))
}
