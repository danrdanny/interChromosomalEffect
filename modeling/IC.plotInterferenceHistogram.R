require(gridExtra)
require(ggplot2)
require(grid)
require(cowplot) # for doing multi-graph plots

data<-read.table("~/projects/GitHub/interChromosomalEffect/modeling/out_bootstrapDCOs.tsv", header=T, sep="\t")

# removed
# geom_histogram(binwidth=100000, colour="black", aes(y=..density.., fill=..count..)) +

# wt is 8344960, IC is 8478101
tmpData <- filter(data,data$chr=="chrX")
chrX.plot <- ggplot(tmpData, aes(x=pos)) +
      stat_function(fun=dnorm,color="red",args=list(mean=mean(tmpData$pos), sd=sd(tmpData$pos))) +
      geom_vline(xintercept=8344960,color="black") +
      geom_vline(xintercept=8478101,color="orange") +
      theme_bw() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

# wt is 11792279, IC is 9811734
tmpData <- filter(data,data$chr=="chr2L")
chr2L.plot <- ggplot(tmpData, aes(x=pos)) +
      stat_function(fun=dnorm,color="red",args=list(mean=mean(tmpData$pos), sd=sd(tmpData$pos))) +
      geom_vline(xintercept=11792279,color="black") +
      geom_vline(xintercept=9811734,color="orange") +
      theme_bw() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

# wt is 10850377, IC is 8448263
tmpData <- filter(data,data$chr=="chr2R")
chr2R.plot <- ggplot(tmpData, aes(x=pos)) +
      stat_function(fun=dnorm,color="red",args=list(mean=mean(tmpData$pos), sd=sd(tmpData$pos))) +
      geom_vline(xintercept=10850377,color="black") +
      geom_vline(xintercept=8448263,color="orange") +
      theme_bw() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

grid.arrange(chr2L.plot, chr2R.plot, chrX.plot, ncol=2)

plot_grid(chrX.plot, chr2L.plot, chr2R.plot, rel_heights=c(1,1), labels=c("chrX","chr2L","chr2R"), nrow=2, ncol=2)
