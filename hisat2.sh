#!/bin/bash
REF=/cluster/home/cxchen/genome/MSU_rice/all.chrs.con.fa
GTF=/cluster/home/cxchen/genome/MSU_rice/all.gtf
GFF=/cluster/home/cxchen/genome/MSU_rice/all.gff3
DIR=/cluster/home/cxchen/Project/HuaweiXu/RNAseq
data=/cluster/sharedata/zhujiankang/Duan/XHW-ASI1-RNAseq-2018-0425/Cleandata
chrom_size=/cluster/home/cxchen/tool/IGVTools/genomes/osativa_7.chrom.sizes
#extract_exons.py $GTF >MSU_rice.exons
#extract_splice_sites.py $GTF >MSU_rice.ss
#hisat2-build -p 28 $REF --ss MSU_rice.ss --exon MSU_rice.exons MSU_rice

##nohup bash /cluster/home/cxchen/Project/HuaweiXu/RNAseq/hisat2.sh OsASi-1 OsASi-2 WT-1 WT-2 2>hisat.log &

mkdir $DIR/FASTQC1
mkdir $DIR/FASTQC2
mkdir $data/tmp

for sample in $@
do
#/cluster/apps/fastqc/0.11.7/fastqc $data/$sample/$sample.R1.fq.gz  $data/$sample/$sample.R2.fq.gz -o $DIR/FASTQC1
#/cluster/apps/TrimGalore/0.4.5/trim_galore --paired $data/$sample/$sample.R1.fq.gz $data/$sample/$sample.R2.fq.gz -o $data/tmp
#/cluster/apps/fastqc/0.11.7/fastqc $data/tmp/$sample.R1_val_1.fq.gz $data/tmp/$sample.R2_val_2.fq.gz -o $DIR/FASTQC2
#hisat2 -p 28 -x MSU_rice -1 $data/tmp/$sample.R1_val_1.fq.gz -2 $data/tmp/$sample.R2_val_2.fq.gz -S $data/tmp/$sample.sam --no-mixed --no-discordant
#/cluster/apps/samtoosl/1.7/bin/samtools view -@ 28 -bS $data/tmp/$sample.sam -o $data/tmp/$sample.bam
#/cluster/apps/samtoosl/1.7/bin/samtools sort -@ 28 $data/tmp/$sample.bam -o $data/tmp/$sample.sorted.bam
#/cluster/apps/samtoosl/1.7/bin/samtools index -@ 28 $data/tmp/$sample.sorted.bam
##infer_experiment.py -i $data/tmp/$sample.sorted.bam -r /cluster/home/cxchen/genome/tair10/TAIR10_GFF3_genes_transposons_noChr.bed
#/cluster/apps/samtoosl/1.7/bin/samtools sort -@ 28  -n $data/tmp/$sample.bam -o $data/tmp/$sample.name.bam
#/cluster/apps/samtoosl/1.7/bin/samtools view -@ 28 -h $data/tmp/$sample.name.bam -o $data/tmp/$sample.name.sam

#igvtools count $data/tmp/$sample.sorted.bam -z 5 -w 25 $DIR/$sample.tdf $chrom_size

htseq-count -f sam -r name -q $data/tmp/$sample.name.sam $GTF > $data/tmp/$sample.htseq_out.txt
########################################################################
samtools depth $data/tmp/$sample.sorted.bam -Q 10 > $DIR/$sample.depth
TOTAL=`cut -f 3 $DIR/$sample.depth |echo $[ $(tr '\n' '+' ) 0 ]`
TA=64107641
########################################################################
sed -n "1,55985p" $data/tmp/$sample.htseq_out.txt >$DIR/$sample.htseq_out.R
 
########################################################################
echo 'all work finished here'
done

