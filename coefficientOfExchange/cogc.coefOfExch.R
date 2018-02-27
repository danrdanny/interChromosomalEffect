require(plyr)
require(dplyr)
require(gridExtra) # for doing multi-graph plots
require(ggplot2)

data<-read.table("~/projects/interChrom/coefficientOfExchange/coefficentOfExchange.coOnly.ice.dat", header=T, sep="\t")
data$Year<-"InterChrom"
data2016<-read.table("projects/interChrom/coefficientOfExchange/coefficentOfExchange.coOnly.wt.dat", header=T, sep="\t")
data2016$Year<-"WT"


data <- rbind(data,data2016)

chrXData <- filter(data,data$Chr=="chrX")
chr2LData <- filter(data,data$Chr=="chr2L")
chr2RData <- filter(data,data$Chr=="chr2R")
#chr3LData <- filter(data,data$Chr=="chr3L")
#chr3RData <- filter(data,data$Chr=="chr3R")

# 2016 and ice data
chrXplot <- ggplot(chrXData, aes(x = Position, y = COEF, color = Year)) + geom_point(size=1) + geom_smooth(alpha=0.1) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(-.01, .2)) + ggtitle("chrX") + ylab("Coefficent of Exchange") + xlab("Pos (Mb)")
chr2Lplot <- ggplot(chr2LData, aes(x=Position, y=COEF, color = Year)) + geom_point(size=1) + geom_smooth(alpha=0.1) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(-.01, .2)) + ggtitle("chr2L") + ylab("Coefficent of Exchange") + xlab("Pos (Mb)")
chr2Rplot <- ggplot(chr2RData, aes(x=Position, y=COEF, color = Year)) + geom_point(size=1) + geom_smooth(alpha=0.1) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(-.01, .2)) + ggtitle("chr2R") + ylab("Coefficent of Exchange") + xlab("Pos (Mb)")
grid.arrange(chr2Lplot, chr2Rplot, chrXplot, ncol=2)


## cm/mb ice data
data<-read.table("~/projects/interChrom/coefficientOfExchange/out_cmMbintervals.ice.tsv", header=T, sep="\t")
data$Year<-"InterChrom"
data2016<-read.table("~/projects/interChrom/coefficientOfExchange/out_cmMbintervals.wt.tsv", header=T, sep="\t")
data2016$Year<-"WT"
data <- rbind(data,data2016)
chrXData <- filter(data,data$Chr=="chrX")
chr2LData <- filter(data,data$Chr=="chr2L")
chr2RData <- filter(data,data$Chr=="chr2R")
chrXplot <- ggplot(chrXData, aes(x = Position, y = COEF, color = Year)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 15)) + ggtitle("chrX") + ylab("cM/Mb") + xlab("Pos (Mb)")
chr2Lplot <- ggplot(chr2LData, aes(x=Position, y=COEF, color = Year)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 15)) + ggtitle("chr2L") + ylab("cM/Mb") + xlab("Pos (Mb)")
chr2Rplot <- ggplot(chr2RData, aes(x=Position, y=COEF, color = Year)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 15)) + ggtitle("chr2R") + ylab("cM/Mb") + xlab("Pos (Mb)")
grid.arrange(chr2Lplot, chr2Rplot, chrXplot, ncol=2)


## cm/mb NCOGC - recombining chromosomes
data<-read.table("~/projects/GitHub/interChromosomalEffect/coefficientOfExchange/out_cmMbintervals.NCOGC.recombining.tsv", header=T, sep="\t")
chrXData <- filter(data,data$Chr=="chrX")
chr2LData <- filter(data,data$Chr=="chr2L")
chr2RData <- filter(data,data$Chr=="chr2R")
chrXplot  <- ggplot(chrXData,  aes(x=Position, y=COEF, color = type)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 6)) + ggtitle("chrX") + ylab("cM/Mb") + xlab("Pos (Mb)")
chr2Lplot <- ggplot(chr2LData, aes(x=Position, y=COEF, color = type)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 6)) + ggtitle("chr2L") + ylab("cM/Mb") + xlab("Pos (Mb)")
chr2Rplot <- ggplot(chr2RData, aes(x=Position, y=COEF, color = type)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 6)) + ggtitle("chr2R") + ylab("cM/Mb") + xlab("Pos (Mb)")
grid.arrange(chr2Lplot, chr2Rplot, chrXplot, ncol=2)

## cm/mb NCOGC - balancer chromosomes
data<-read.table("~/projects/GitHub/interChromosomalEffect/coefficientOfExchange/out_cmMbintervals.NCOGC.balancers.tsv", header=T, sep="\t")
chr2LData <- filter(data,data$Chr=="chr2L")
chr2RData <- filter(data,data$Chr=="chr2R")
chr3LData <- filter(data,data$Chr=="chr3L")
chr3RData <- filter(data,data$Chr=="chr3R")

chr2Lplot <- ggplot(chr2LData, aes(x=Position, y=COEF, color = type)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 10)) + ggtitle("chr2L") + ylab("cM/Mb") + xlab("Pos (Mb)")
chr2Rplot <- ggplot(chr2RData, aes(x=Position, y=COEF, color = type)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 10)) + ggtitle("chr2R") + ylab("cM/Mb") + xlab("Pos (Mb)")
chr3Lplot <- ggplot(chr3LData, aes(x=Position, y=COEF, color = type)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 10)) + ggtitle("chr3L") + ylab("cM/Mb") + xlab("Pos (Mb)")
chr3Rplot <- ggplot(chr3RData, aes(x=Position, y=COEF, color = type)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 10)) + ggtitle("chr3R") + ylab("cM/Mb") + xlab("Pos (Mb)")
grid.arrange(chr2Lplot, chr2Rplot, chr3Lplot, chr3Rplot, ncol=2)

## cm/mb NCOGC onto balancer chromosomes after rearranging them
data<-read.table("~/projects/GitHub/interChromosomalEffect/coefficientOfExchange/out_cmMbintervals.NCOGC.balancerOrder.tsv", header=T, sep="\t")
chr2LData <- filter(data,data$Chr=="chr2L")
chr2RData <- filter(data,data$Chr=="chr2R")
chr3LData <- filter(data,data$Chr=="chr3L")
chr3RData <- filter(data,data$Chr=="chr3R")

chr2Lplot <- ggplot(chr2LData, aes(x=Position, y=COEF, color = type)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 10)) + ggtitle("chr2L") + ylab("cM/Mb") + xlab("Pos (Mb)")
chr2Rplot <- ggplot(chr2RData, aes(x=Position, y=COEF, color = type)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 10)) + ggtitle("chr2R") + ylab("cM/Mb") + xlab("Pos (Mb)")
chr3Lplot <- ggplot(chr3LData, aes(x=Position, y=COEF, color = type)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,5000000,10000000,15000000,20000000,25000000,30000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 10)) + ggtitle("chr3L") + ylab("cM/Mb") + xlab("Pos (Mb)")
chr3Rplot <- ggplot(chr3RData, aes(x=Position, y=COEF, color = type)) + geom_point(size=1) + geom_smooth(alpha=0) + theme_minimal() + scale_x_continuous(breaks = c(0,10000000,20000000,30000000,40000000)) + scale_y_continuous(expand = c(0, 0), limit = c(0, 10)) + ggtitle("chr3R") + ylab("cM/Mb") + xlab("Pos (Mb)")
grid.arrange(chr2Lplot, chr2Rplot, chr3Lplot, chr3Rplot, ncol=2)




## Unused below here


# 2016 plots
chrXplot <- ggplot(chrXData, aes(x=Position, y=COEF)) + geom_point(size=3) + geom_smooth(alpha=1) + theme_minimal() + scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) + scale_y_continuous(expand = c(0, 0), limit = c(-.05, .25)) + ggtitle("chrX") + ylab("Coefficent of Exchange") + xlab("Pos (Mb)") + theme(legend.position="none")
chr2Lplot <- ggplot(chr2LData, aes(x=Position, y=COEF)) + geom_point(size=3) + geom_smooth(alpha=1) + theme_minimal() + scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) + scale_y_continuous(expand = c(0, 0), limit = c(-.05, .25)) + ggtitle("chr2L") + ylab("Coefficent of Exchange") + xlab("Pos (Mb)") + theme(legend.position="none")
chr2Rplot <- ggplot(chr2RData, aes(x=Position, y=COEF)) + geom_point(size=3) + geom_smooth(alpha=1) + theme_minimal() + scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) + scale_y_continuous(expand = c(0, 0), limit = c(-.05, .25)) + ggtitle("chr2R") + ylab("Coefficent of Exchange") + xlab("Pos (Mb)") + theme(legend.position="none")
chr3Lplot <- ggplot(chr3LData, aes(x=Position, y=COEF)) + geom_point(size=3) + geom_smooth(alpha=1) + theme_minimal() + scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) + scale_y_continuous(expand = c(0, 0), limit = c(-.05, .25)) + ggtitle("chr3L") + ylab("Coefficent of Exchange") + xlab("Pos (Mb)") + theme(legend.position="none")
chr3Rplot <- ggplot(chr3RData, aes(x=Position, y=COEF)) + geom_point(size=3) + geom_smooth(alpha=1) + theme_minimal() + scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) + scale_y_continuous(expand = c(0, 0), limit = c(-.05, .25)) + ggtitle("chr3R") + ylab("Coefficent of Exchange") + xlab("Pos (Mb)") + theme(legend.position="none")

grid.arrange(chr2Lplot, chr2Rplot, chrXplot, ncol=2)

## Below is unused, but can be adapted to plot COs and NCOs together.

data<-read.table("~/Dropbox/CO and GC Paper/data/coefficentOfExchange.coandgc.dat", header=T, sep="\t")
data<-read.table("~/Dropbox/CO and GC Paper/data/coefficentOfExchange.all.dat", header=T, sep="\t")
data<-read.table("~/Dropbox/CO and GC Paper/data/coefficentOfExchange.cogcsep.dat", header=T, sep="\t")

chrXData <- filter(data,data$Chr=="X")
chr2LData <- filter(data,data$Chr=="2L")
chr2RData <- filter(data,data$Chr=="2R")
chr3LData <- filter(data,data$Chr=="3L")
chr3RData <- filter(data,data$Chr=="3R")

chrXplot <- ggplot(chrXData, aes(x=Band, y=coef, color=cond)) + geom_point(aes(shape=factor(cond)),size=3) + geom_smooth() + theme_minimal() + scale_y_continuous(expand = c(0, 0), limit = c(-.1, .3)) + ggtitle("chrX") + ylab("Coefficent of Exchange") + theme(legend.position="none")
chr2Lplot <- ggplot(chr2LData, aes(x=Band, y=coef, color=cond)) + geom_point(aes(shape=factor(cond)),size=3) + geom_smooth() + theme_minimal() + scale_y_continuous(expand = c(0, 0), limit = c(-.1, .3)) + ggtitle("chr2L") + ylab("Coefficent of Exchange") + theme(legend.position="none")
chr2Rplot <- ggplot(chr2RData, aes(x=Band, y=coef, color=cond)) + geom_point(aes(shape=factor(cond)),size=3) + geom_smooth() + theme_minimal() + scale_y_continuous(expand = c(0, 0), limit = c(-.1, .3)) + ggtitle("chr2R") + ylab("Coefficent of Exchange") + theme(legend.position="none")
chr3Lplot <- ggplot(chr3LData, aes(x=Band, y=coef, color=cond)) + geom_point(aes(shape=factor(cond)),size=3) + geom_smooth() + theme_minimal() + scale_y_continuous(expand = c(0, 0), limit = c(-.1, .3)) + ggtitle("chr3L") + ylab("Coefficent of Exchange") + theme(legend.position="none")
chr3Rplot <- ggplot(chr3RData, aes(x=Band, y=coef, color=cond)) + geom_point(aes(shape=factor(cond)),size=3) + geom_smooth() + theme_minimal() + scale_y_continuous(expand = c(0, 0), limit = c(-.1, .3)) + ggtitle("chr3R") + ylab("Coefficent of Exchange") + theme(legend.position="none")

chrXplot <- ggplot(chrXData, aes(x=Band, y=coef, color=cond)) + geom_point(aes(shape=factor(cond)),size=3) + geom_smooth() + theme_minimal() + ggtitle("chrX") + ylab("Coefficent of Exchange") + theme(legend.position="none")
chr2Lplot <- ggplot(chr2LData, aes(x=Band, y=coef, color=cond)) + geom_point(aes(shape=factor(cond)),size=3) + geom_smooth() + theme_minimal() + ggtitle("chr2L") + ylab("Coefficent of Exchange") + theme(legend.position="none")
chr2Rplot <- ggplot(chr2RData, aes(x=Band, y=coef, color=cond)) + geom_point(aes(shape=factor(cond)),size=3) + geom_smooth() + theme_minimal() + ggtitle("chr2R") + ylab("Coefficent of Exchange") + theme(legend.position="none")
chr3Lplot <- ggplot(chr3LData, aes(x=Band, y=coef, color=cond)) + geom_point(aes(shape=factor(cond)),size=3) + geom_smooth() + theme_minimal() + ggtitle("chr3L") + ylab("Coefficent of Exchange") + theme(legend.position="none")
chr3Rplot <- ggplot(chr3RData, aes(x=Band, y=coef, color=cond)) + geom_point(aes(shape=factor(cond)),size=3) + geom_smooth() + theme_minimal() + ggtitle("chr3R") + ylab("Coefficent of Exchange") + theme(legend.position="none")

