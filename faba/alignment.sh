
 STAR --runThreadN 16 --genomeDir STAR_utr --readFilesIn ../20251230_atlas_3lanes/BMK250917-CR226-ZX01-0101/BMK_DATA_20251227193551_1/merged_data/36HAI-1_R2_001.fastq.gz ../20251230_atlas_3lanes/BMK250917-CR226-ZX01-0101/BMK_DATA_20251227193551_1/merged_data/36HAI-1_R1_001.fastq.gz \
    --readFilesCommand zcat --soloType CB_UMI_Simple --soloCBwhitelist 3M-3pgex-may-2023.txt --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIlen 12 --soloBarcodeReadLength 0\
    --soloCBlen 16 --soloCBstart 1 --soloUMIstart 17  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix star_faba_intron_r2_36hair1/ --outFilterMismatchNmax 2 \
    --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --soloFeatures GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS\
    --soloMultiMappers EM --outReadsUnmapped Fastx --soloOutFileNames star_faba_intron_r2_36hair1/ features.tsv barcodes.tsv matrix.mtx 

#change matrix.mtx name to get the matrix including multimap
#empty cell filter
#STAR --runMode soloCellFiltering  star_faba_iseq_intron/star_faba_iseq_intron/GeneFull/raw   star_faba_iseq_intron/star_faba_iseq_intron/GeneFull/filter_cr/   --soloCellFilter EmptyDrops_CR

 
