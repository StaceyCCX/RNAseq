rm (list=ls())

setwd("/cluster/home/cxchen/Project/HuaweiXu/RNAseq")

wt1=read.table("WT-1.htseq_out_TE_up_down.R")
wt2=read.table("WT-2.htseq_out_TE_up_down.R")

OsASi_1=read.table("OsASi-1.htseq_out_TE_up_down.R")
OsASi_2=read.table("OsASi-2.htseq_out_TE_up_down.R")

TE_up_down=data.frame(wt1=wt1$V2,wt2=wt2$V2,OsASi_1=OsASi_1$V2,OsASi_2=OsASi_2$V2)
rownames(TE_up_down)=wt1$V1

group=c(rep("wt",2),rep("OsASi",2))
library(edgeR)
cds = DGEList(TE_up_down,group =group)
cds = calcNormFactors(cds)
keep = rowSums(cpm(cds)>1)>=2
cds =cds[keep,]
log2counts = log2(cds$counts+1)
pdf ("TE_up_down_reads_wt_OsASi_correlation.pdf")
panel.cor <- function(x,y, ...)
{ 
par(usr=c(0,1,0,1))
txt <- as.character(format(cor(x,y,method=c("spearman")),digits=4))
text(0.5,0.5,txt,cex =2*abs(cor(x,y)))
}
pairs(log2counts[,1:4],upper.panel=panel.cor,,main="Relationship of wt replicates and OsASi replicates in Oryza sativa")
dev.off()

cds=estimateCommonDisp(cds,verbose=TRUE)
cds = estimateTagwiseDisp(cds)
cpmcounts = cpm(cds$counts+1)
tcounts = t(cpmcounts)
pca.total = prcomp(log2(tcounts), retx=TRUE)
pdf("TE_up_down_reads_wt_OsASi_pca.pdf")
c=round(100*summary(pca.total)$importance[2,1],digits=2)
d=round(100*summary(pca.total)$importance[2,2],digits=2)
tt=c("wt1","wt2","OsASi_1","OsASi_2")
plot(pca.total$x[,1:2], pch=c(15:17), xlab=paste("PC1(",c,"% Proportion of Variance)"),ylab=paste("PC2(",d,"%) Proportion of Variance"),col=c(rep("black",2),rep("red",2)),main="PCA Plot of Samples")
legend("bottomright",cex=0.6,border=F, legend=tt,pch=c(15:17), col=c(rep("black",2),rep("red",2)),bty="n")
dev.off()
de.out=exactTest(cds,pair=c("wt","OsASi"))
rawcount = as.matrix(TE_up_down)[(rownames(topTags(de.out,n=100000))),]
rawcount=rawcount[order(rownames(rawcount)),]
table=topTags(de.out,n=100000)$table
table=table[order(rownames(table)),]
sam_cpm=cpm(cds$count)[row.names(rawcount),]
TE_up_down=cbind(rawcount,sam_cpm,table)

TE_up_down=TE_up_down[order(TE_up_down$PValue),]
write.table(TE_up_down,"TE_up_down_reads_wt_and_OsASi.edgeR.out",sep="\t",row.names=T,col.names=F,quote=F)

de.up=subset(TE_up_down,((TE_up_down$FDR<0.05) & (TE_up_down$logFC) >=1))
de.down=subset(TE_up_down,((TE_up_down$FDR<0.05) & (TE_up_down$logFC) <= -1))
write.table(de.up,"TE_up_down_reads_ee_DEG.up",sep="\t",row.names=T,col.names=F,quote=F)
write.table(de.down,"TE_up_down_reads_ee_DEG.down",sep="\t",row.names=T,col.names=F,quote=F)

c1up=nrow(subset(TE_up_down,((TE_up_down$FDR<0.05) & (TE_up_down$logFC) >=1)))
c1down=nrow(subset(TE_up_down,((TE_up_down$FDR<0.05) & (TE_up_down$logFC) <= -1)))
c0.58up=nrow(subset(TE_up_down,((TE_up_down$FDR<0.05) & (TE_up_down$logFC) >=0.58)))
c0.58down=nrow(subset(TE_up_down,((TE_up_down$FDR<0.05) & (TE_up_down$logFC) <= -0.58)))
de.log=matrix(c(c1up,c1down,c0.58up,c0.58down),nrow=2)
write.table(de.log,"TE_up_down_reads_ee_Cutoff.log",sep="\t",row.names=F,col.names=F,quote=F)
