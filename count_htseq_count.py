import sys 
import string

with open('/cluster/home/cxchen/genome/MSU_rice/all.gff3') as all_gff3:
	id=[]
	strand=[]
	gene={}
	for line in all_gff3:
		if len(line.split('\t'))<3:
			pass
		else:
			if line.split('\t')[2].find('gene')>=0:
				id=line.split('\t')[8].split(';')[0].strip('ID=')
				strand=line.split('\t')[6]
				gene[id]=strand
with open(sys.argv[1]) as file_1:
	id=[]
	reads=[]
	gene_1={}
	for line in file_1:
		id=line.split('\t')[0]
		reads=line.split('\t')[1].strip('\n')
		gene_1[id]=reads
with open(sys.argv[2]) as file_2:
	id=[]
	reads=[]
	gene_2={}
	for line in file_2:
		id=line.split('\t')[0]
		reads=line.split('\t')[1].strip('\n')
		gene_2[id]=reads
with open(sys.argv[3],'w') as result:
	for item in gene:
		item1=item+'.1'
		item2=item+'.2'
		if item1 in gene_1 and item2 in gene_2:
			if gene[item]=='+'  and float(gene_1.get(item1))>0:
				result.write(item+'\t'+gene_1.get(item1)+'\t'+gene_2.get(item2)+'\t'+str(float(gene_2.get(item2))/float(gene_1.get(item1)))+'\n')
			if gene[item]=='-' and float(gene_2.get(item2))>0:
				result.write(item+'\t'+gene_1.get(item1)+'\t'+gene_2.get(item2)+'\t'+str(float(gene_1.get(item1))/float(gene_2.get(item2)))+'\n')
#python /cluster/home/cxchen/Project/HuaweiXu/RNAseq/count_htseq_count.py OsASi-1_htseq_count_TE_cantaing_gene.1.txt OsASi-1_htseq_count_TE_cantaing_gene.2.txt OsASi-1_htseq_count_TE_cantaing_gene.txt
#this script use to calculate intronic TE downstream/upstream ratio from htseq count out result
#note:
	#we should set cutoff so that both reads in 5' and 3' not less than 20
	#but for treatment,just keep 5'>0 (divisor must not equal to 0)
