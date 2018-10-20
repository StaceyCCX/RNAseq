#PBS -N hisat_01
#PBS -q batch
#PBS -l nodes=30:ppn=20
#PBS -j oe
#PBS -t 1-100

# use Differncial expression RNA seq method to idetified RIP peaks

dir=/cluster/sharedata/zhujiankang/Duan/ASI1_RIP

ARRAY=(ASI1-RIP-1  ASI1-RIP-2  ASI1-RIP-input  WT-RIP-1  WT-RIP-2  WT-RIP-input)
mkdir $dir/BAM/


for sample in $@
do

sickle pe -f ${ARRAY[$PBS_ARRAYID]}/*L001_R1_001.fastq.gz -r ${ARRAY[$PBS_ARRAYID]}/*L001_R2_001.fastq.gz  -t sanger -l 50 -q 20 -o $dir/trimed_${ARRAY[$PBS_ARRAYID]}.L001_R1.fastq.gz -p $dir/trimed_${ARRAY[$PBS_ARRAYID]}.L001_R2.fastq.gz -s $dir/trimed_${ARRAY[$PBS_ARRAYID]}.L001_single.fastq.gz &
sickle pe -f ${ARRAY[$PBS_ARRAYID]}/*L002_R1_001.fastq.gz -r ${ARRAY[$PBS_ARRAYID]}/*L002_R2_001.fastq.gz  -t sanger -l 50 -q 20 -o $dir/trimed_${ARRAY[$PBS_ARRAYID]}.L002_R1.fastq.gz -p $dir/trimed_${ARRAY[$PBS_ARRAYID]}.L002_R2.fastq.gz -s $dir/trimed_${ARRAY[$PBS_ARRAYID]}.L002_single.fastq.gz &

######################## build genome index ##############################
STAR --runMode genomeGenerate --runThreadN 20 --genomeFastaFiles Chr_TAIR10_chr_all.fas --sjdbGTFfile Araport11_GFF3_genes_transposons.201606.gff2.gtf --genomeDir STARindex_transcripts/

STAR --runThreadN 10 --genomeDir /cluster/home/cxchen/genome/tair10/STARindex_transcripts/ --readFilesIn $dir/fastq/trimed_${ARRAY[$PBS_ARRAYID]}.R1.fastq.gz $dir/fastq/trimed_${ARRAY[$PBS_ARRAYID]}.R2.fastq.gz --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $dir/BAM/${ARRAY[$PBS_ARRAYID]}.  -outWigType wiggle

gffread -w Araport11_GFF3_genes_transposons.201606.gff2.gtf.fas -g Chr_TAIR10_chr_all.fas Araport11_GFF3_genes_transposons.201606.gff2.gtf
salmon quant -t /cluster/home/cxchen/genome/tair10/Araport11_GFF3_genes_transposons.201606.gff2.gtf.fas -l A -a $dir/BAM/${ARRAY[$PBS_ARRAYID]}.Aligned.toTranscriptome.out.bam -o ${ARRAY[$PBS_ARRAYID]}.quant

done
