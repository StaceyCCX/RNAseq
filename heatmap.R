rm(list=ls())

setwd("~/myProject/other/XSS")
library(gplots)
library(d3heatmap)
library(ComplexHeatmap)

data=read.table('abiotic_lixt.txt',header = F)
df=data[,2:5]
row=data[,1]
column=c('wt','wt','28-1','28-1')
dimnames(df)=list(row,column)
df=as.matrix(df)
par(mar=c(1,0.5,1,2))
heatmap.2(df, scale = "row", col=bluered(100), trace = "none", density.info = "none")

