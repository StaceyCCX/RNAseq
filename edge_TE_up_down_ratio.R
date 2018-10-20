rm (list=ls())

setwd("/cluster/home/cxchen/Project/HuaweiXu/RNAseq")

wt1=read.table("WT-1_htseq_count_TE_cantaing_gene.R")
wt2=read.table("WT-2_htseq_count_TE_cantaing_gene.R")

asi_1=read.table("OsASi-1_htseq_count_TE_cantaing_gene.R")
asi_2=read.table("OsASi-2_htseq_count_TE_cantaing_gene.R")

all=data.frame(wt1=wt1$V2,wt2=wt2$V2,asi_1=asi_1$V2,asi_2=asi_2$V2)
rownames(all)=(wt1$V1)

group=c(rep("wt",2),rep("asi",2))
library(edgeR)
cds = DGEList(all,group =group)
cds = calcNormFactors(cds)
cds=estimateCommonDisp(cds,verbose=TRUE)
cds = estimateTagwiseDisp(cds)
de.out=exactTest(cds,pair=c("wt","asi"))
rawcount = as.matrix(all)[(rownames(topTags(de.out,n=100000))),]
rawcount=rawcount[order(rownames(rawcount)),]
table=topTags(de.out,n=100000)$table
table=table[order(rownames(table)),]
sam_cpm=cpm(cds$count)[row.names(rawcount),]
all=cbind(rawcount,sam_cpm,table)

all=all[order(all$PValue),]
write.table(all,"TE_cantaing_gene_0.58_wt_and_asi_ratio.edgeR.out",sep="\t",row.names=T,col.names=F,quote=F)

de.up=subset(all,((all$logFC) >=0.58))
de.down=subset(all,((all$logFC) <=-0.58))
write.table(de.up,"TE_cantaing_gene_0.58_wt_and_asi_DEG_ratio.up",sep="\t",row.names=T,col.names=F,quote=F)
write.table(de.down,"TE_cantaing_gene_0.58_wt_and_asi_DEG_ratio.down",sep="\t",row.names=T,col.names=F,quote=F)


