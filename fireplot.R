rm(list=ls())

setwd("~/myProject/other/XSS")
library(ggthemes)
library(Cairo)
library(ggplot2)

data <- read.delim("all_wt_and_28-1.edgeR.out",header = F,sep="\t",na.strings = "")
data$threshold <- as.factor(ifelse(data$V12 < 0.05 & abs(data$V10) >=1.5,ifelse(data$V10 > 1.5 ,'Up','Down'),'Not')) #V12:pvalues;v10:log2(foldchange)
Cairo(file="all_wt_and_28-1.edgeR.out.png", 
	type="png",
	units="in",
	bg="white",
	width=5.5, 
	height=5, 
	pointsize=12, 
	dpi=300)
ggplot(data=data, 
	aes(x=V10, y =-log10(V12),  
	colour=threshold,fill=threshold)) +
	scale_color_manual(values=c("blue", "grey","red"))+
	geom_point(alpha=0.4, size=1.2) +
	xlim(c(-12, 12)) +
	theme_bw(base_size = 12, base_family = "Times") +
	geom_vline(xintercept=c(-1.5,1.5),lty=4,col="grey",lwd=0.6)+
	geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
	theme(legend.position="right",
	panel.grid=element_blank(),
	legend.title = element_blank(),
	legend.text= element_text(face="bold", color="black",family = "Times", size=8),
	plot.title = element_text(hjust = 0.5),
	axis.text.x = element_text(face="bold", color="black", size=12),
	axis.text.y = element_text(face="bold",  color="black", size=12),
	axis.title.x = element_text(face="bold", color="black", size=12),
	axis.title.y = element_text(face="bold",color="black", size=12))+
	labs(x="log2 (fold change)",y="-log10 (p-value)",title="Volcano picture of DEG")
dev.off()
