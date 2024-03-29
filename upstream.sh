#prefetch SRR23341322 
#fastq-dump --split-file --gzip ./SRR23341322/SRR23341322.sra -O ./fastq/
#prepare for cell ranger
cat SRR_Acc_List-9245-3.txt | while read i ;do (mv ${i}_1*.gz ${i}_S1_L001_I1_001.fastq.gz;mv ${i}_2*.gz ${i}_S1_L001_I2_001.fastq.gz;mv ${i}_3*.gz ${i}_S1_L001_R1_001.fastq.gz;mv ${i}_4*.gz ${i}_S1_L001_R2_001.fastq.gz);done

#create ref https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr
# check the number of gene_biotype
cat MtrunA17r5.0-ANR-EGN-r1.9.gtf |grep -v "#" |awk -v FS='gene_biotype ' 'NF>1{print $2}'|awk -F ";" '{print $1}'|sort | uniq -c
#   1214 "miRNA"
# 524306 "mRNA"
#  18469 "ncRNA"
#   2397 "pre_miRNA"
# 735219 "repeat_region"
#    186 "rRNA"
#   2922 "tRNA"
##filter miRNA because it have no exon
cellranger mkgtf \
  MtrunA17r5.0-ANR-EGN-r1.9.gtf \
   MtrunA17r5.0-ANR-EGN-r1.9.filtered.gtf \
  --attribute=gene_biotype:mRNA \
  --attribute=gene_biotype:ncRNA \
  --attribute=gene_biotype:pre_miRNA \
  --attribute=gene_biotype:repeat_region \
  --attribute=gene_biotype:rRNA \
  --attribute=gene_biotype:tRNA
cellranger mkref   --genome=Mtrun   --fasta=/home/mlei/ljscseq/medicago_truncatula/genome/MtrunA17r5.0-20161119-ANR.genome.fasta   --genes=/home/mlei/ljscseq/medicago_truncatula/genome/MtrunA17r5.0-ANR-EGN-r1.9.filtered.gtf
cellranger count --id=result_SRR23341322 \
                  --transcriptome=/home/mlei/ljscseq/medicago_truncatula/scrna-seq/Mtrun \
                  --fastqs=/home/mlei/ljscseq/medicago_truncatula/scrna-seq/fastq/ \
                  --sample=SRR23341322




