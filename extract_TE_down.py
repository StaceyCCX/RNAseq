#python extract_TE_down.py TE_up_down_reads_ee_DEG.down result_TE_up_down_reads_ee_DEG.down
#this script use to extract intron-containing genes which shows down stream of intron-TE was down regulated and up stream of intron-TE was not down regulated.
import sys
import string

gene_id=[]
strand=[]
gene={}

with open('/cluster/home/cxchen/genome/MSU_rice/all.gff3','r') as gene_file:
	for line in gene_file:
		if line.startswith('#'):
			pass
		if len(line.split('\t'))<3:
			pass
		else:
			if line.split('\t')[2]=='gene':
				gene_id=line.split('\t')[8].split(';')[0].strip('ID=')
				strand=line.split('\t')[6]
				gene[gene_id]=strand
id=[]
wait_file=open(sys.argv[1]).readlines()
for line in wait_file:
	id.append(line.split('\t')[0])
TE_down=[]
for item in gene:
	if gene[item]=='+':
		if item+'.2' in id and item+'.1' not in id:
			TE_down.append(item)
	if gene[item]=='-':
		if item+'.1' in id and item+'.2' not in id:
			TE_down.append(item)
with open(sys.argv[2],'w') as TE_down_file:
	for line in wait_file:
		if int(line.split('\t')[1])>=20 and int(line.split('\t')[2])>=20:
			if line.split('\t')[0].split('.')[0] in TE_down:
				line=line.replace(line.split('\t')[0],line.split('\t')[0].split('.')[0])
				TE_down_file.write(line)	
