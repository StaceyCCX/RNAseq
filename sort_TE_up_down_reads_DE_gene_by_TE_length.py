import sys
import string

gene_id=[]
TE_length=[]
TE={}
with open('/cluster/home/cxchen/genome/MSU_rice/intron_extract/inTE_wo.bed.sort.result.rmdup_TElength.count.sort','r') as TE_length_file:
	for line in TE_length_file:
		gene_id=line.split('\t')[9]
		TE_length=line.split('\t')[12]
		TE[gene_id]=TE_length

DE_gene_TElength=open(sys.argv[2],'w')
with open(sys.argv[1],'r') as DE_gene_file:
	for line in DE_gene_file:
		id=line.split('\t')[0] 
		DE_gene_TElength.write(line.strip('\n')+'\t'+TE[id])
DE_gene_TElength.close()


			
