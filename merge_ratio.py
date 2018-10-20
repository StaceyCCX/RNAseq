import sys
import string

with open('WT-1_htseq_count_TE_cantaing_gene.txt') as WT_1:
	id=[]
	ratio=[]
	wt1={}
	for line in WT_1:
		id=line.split('\t')[0]
		ratio=line.split('\t')[1]+'\t'+line.split('\t')[2]+'\t'+line.split('\t')[3].strip('\n')
		wt1[id]=ratio

with open('WT-2_htseq_count_TE_cantaing_gene.txt') as WT_2:
	id=[]
	ratio=[]
	wt2={}
	for line in WT_2:
		id=line.split('\t')[0]
		ratio=line.split('\t')[1]+'\t'+line.split('\t')[2]+'\t'+line.split('\t')[3].strip('\n')
		wt2[id]=ratio

with open('OsASi-1_htseq_count_TE_cantaing_gene.txt') as OsASi_1:
	id=[]
	ratio=[]
	osasi1={}
	for line in OsASi_1:
		id=line.split('\t')[0]
		ratio=line.split('\t')[1]+'\t'+line.split('\t')[2]+'\t'+line.split('\t')[3].strip('\n')
		osasi1[id]=ratio

with open('OsASi-2_htseq_count_TE_cantaing_gene.txt') as OsASi_2:
	id=[]
	ratio=[]
	osasi2={}
	for line in OsASi_2:
		id=line.split('\t')[0]
		ratio=line.split('\t')[1]+'\t'+line.split('\t')[2]+'\t'+line.split('\t')[3].strip('\n')
		osasi2[id]=ratio
common=wt1.keys()&wt2.keys()&osasi1.keys()&osasi2.keys()

Wwt1=open('WT-1_htseq_count_TE_cantaing_gene.R','w')
Wwt2=open('WT-2_htseq_count_TE_cantaing_gene.R','w')
Wosasi1=open('OsASi-1_htseq_count_TE_cantaing_gene.R','w')
Wosasi2=open('OsASi-2_htseq_count_TE_cantaing_gene.R','w')
total=open('WT-1_WT-2_OsASi-1_OsASi-2_TE_reads_ratio.txt','w')

for item in common:
	Wwt1.write(item+'\t'+wt1[item].split('\t')[2]+'\n')
	Wwt2.write(item+'\t'+wt2[item].split('\t')[2]+'\n')
	Wosasi1.write(item+'\t'+osasi1[item].split('\t')[2]+'\n')
	Wosasi2.write(item+'\t'+osasi2[item].split('\t')[2]+'\n')
	total.write(item+'\t'+wt1[item]+'\t'+wt2[item]+'\t'+osasi1[item]+'\t'+osasi2[item]+'\n')

Wwt1.close()
Wwt2.close()
Wosasi1.close()
Wosasi2.close()
total.close()
