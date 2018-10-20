#PBS -N hisat_01
#PBS -q batch
#PBS -l nodes=5:ppn=28
#PBS -j oe
#PBS -t 1-100

REF=/cluster/home/cxchen/genome/MSU_rice/all.chrs.con.fa
GTF=/cluster/home/cxchen/genome/MSU_rice/all.gtf
modified_GTF=/cluster/oldhome/cxchen/MSU_rice/intron_extract/target_up_down_gene1_rmdup.gtf
DIR=/cluster/home/cxchen/Project/HuaweiXu/RNAseq
data=/cluster/sharedata/zhujiankang/Duan/XHW-ASI1-RNAseq-2018-0425/Cleandata
chrom_size=/cluster/home/cxchen/tool/IGVTools/genomes/osativa_7.chrom.sizes

################### building genome file index ###################################
extract_exons.py GTF >MSU_rice.exons
extract_splice_sites.py GTF >MSU_rice.ss
hisat2-build -p 28 $REF --ss MSU_rice.ss --exon MSU_rice.exons MSU_rice

################## start working #################################################
ARRAY=(OsASi-1 OsASi-2 WT-1 WT-2)

##### data clean and fastqc #####
/cluster/apps/fastqc/0.11.7/fastqc/$data/${ARRAY[$PBS_ARRAYID]}/${ARRAY[$PBS_ARRAYID]}.R1.fq.gz  $data/${ARRAY[$PBS_ARRAYID]}/${ARRAY[$PBS_ARRAYID]}.R2.fq.gz -o $DIR/FASTQC1
/cluster/apps/TrimGalore/0.4.5/trim_galore --paired $data/${ARRAY[$PBS_ARRAYID]}/${ARRAY[$PBS_ARRAYID]}.R1.fq.gz $data/${ARRAY[$PBS_ARRAYID]}/${ARRAY[$PBS_ARRAYID]}.R2.fq.gz -o $data/tmp
/cluster/apps/fastqc/0.11.7/fastqc $data/tmp/${ARRAY[$PBS_ARRAYID]}.R1_val_1.fq.gz $data/tmp/${ARRAY[$PBS_ARRAYID]}.R2_val_2.fq.gz -o $DIR/FASTQC2

##### mapping to genome #####
hisat2 -p 28 -x MSU_rice -1 $data/tmp/${ARRAY[$PBS_ARRAYID]}.R1_val_1.fq.gz -2 $data/tmp/${ARRAY[$PBS_ARRAYID]}.R2_val_2.fq.gz -S $data/tmp/${ARRAY[$PBS_ARRAYID]}.sam --no-mixed --no-discordant
/cluster/apps/samtoosl/1.7/bin/samtools view -@ 28 -bS $data/tmp/${ARRAY[$PBS_ARRAYID]}.sam -q 30 -o $data/tmp/${ARRAY[$PBS_ARRAYID]}.bam
/cluster/apps/samtoosl/1.7/bin/samtools sort -@ 28 $data/tmp/${ARRAY[$PBS_ARRAYID]}.bam -o $data/tmp/${ARRAY[$PBS_ARRAYID]}.sorted.bam
/cluster/apps/samtoosl/1.7/bin/samtools index -@ 28 $data/tmp/${ARRAY[$PBS_ARRAYID]}.sorted.bam

##infer_experiment.py -i $data/tmp/${ARRAY[$PBS_ARRAYID]}.sorted.bam -r /cluster/home/cxchen/genome/tair10/TAIR10_GFF3_genes_transposons_noChr.bed

##### reads count for normal DE genes and modified DE genes #####
/cluster/apps/samtoosl/1.7/bin/samtools sort -@ 28  -n $data/tmp/${ARRAY[$PBS_ARRAYID]}.bam -o $data/tmp/${ARRAY[$PBS_ARRAYID]}.name.bam
/cluster/apps/samtoosl/1.7/bin/samtools view -@ 28 -h $data/tmp/${ARRAY[$PBS_ARRAYID]}.name.bam -o $data/tmp/${ARRAY[$PBS_ARRAYID]}.name.sam
igvtools count $data/tmp/${ARRAY[$PBS_ARRAYID]}.sorted.bam -z 5 -w 25 $DIR/${ARRAY[$PBS_ARRAYID]}.tdf $chrom_size

htseq-count -f sam -r name -q $data/tmp/${ARRAY[$PBS_ARRAYID]}.name.sam $GTF > $data/tmp/htseq_out_${ARRAY[$PBS_ARRAYID]}.txt
sed -n "1,55985p" $data/tmp/htseq_out_${ARRAY[$PBS_ARRAYID]}.txt >$DIR/htseq_out_${ARRAY[$PBS_ARRAYID]}.R

htseq-count -f sam -r name -q $data/tmp/${ARRAY[$PBS_ARRAYID]}.name.sam $modified_GTF > $data/tmp/htseq_out_${ARRAY[$PBS_ARRAYID]}.TE_up_down.txt
sed -n "1,34906p" $data/tmp/htseq_out_${ARRAY[$PBS_ARRAYID]}.TE_up_down.txt >$DIR/htseq_out_${ARRAY[$PBS_ARRAYID]}.TE_up_down.R
 
################### endding ##########################################################
echo 'all work finished'

