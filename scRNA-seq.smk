configfile: "/home/mlei/ljscseq/medicago_truncatula/scRNA-seq.yaml"

import pandas as pd

sra = pd.read_csv("/home/mlei/ljscseq/medicago_truncatula/scrna-seq/srr.txt",sep="\t")
SAMPLES = sra["Run"].tolist()#[1:5]


WORKING_DIR = config['wd']
RESULT_DIR = config["rd"]
RESULT_DIR2 = config["rd"]

REF_dir = config["refdir"]
genome = config["ref"]


rule all:
    input:    
        i1 = expand(WORKING_DIR + "fastq/{sample}_S1_L001_I1_001.fastq.gz",sample = SAMPLES),
        i2= expand(WORKING_DIR + "fastq/{sample}_S1_L001_I2_001.fastq.gz",sample = SAMPLES),
        r1= expand(WORKING_DIR + "fastq/{sample}_S1_L001_R1_001.fastq.gz",sample = SAMPLES),
        r2= expand(WORKING_DIR + "fastq/{sample}_S1_L001_R2_001.fastq.gz",sample = SAMPLES),
        results = expand(WORKING_DIR+ "result_{sample}",sample = SAMPLES)

    
rule get_SRR_files:
    output:
        i1 = WORKING_DIR + "fastq/{sample}_S1_L001_I1_001.fastq.gz",
        i2= WORKING_DIR + "fastq/{sample}_S1_L001_I2_001.fastq.gz",
        r1= WORKING_DIR + "fastq/{sample}_S1_L001_R1_001.fastq.gz",
        r2= WORKING_DIR + "fastq/{sample}_S1_L001_R2_001.fastq.gz"
    params:
       sample = "{sample}",
       DIR = WORKING_DIR+"fastq/"
#       tmpd = WORKING_DIR + "tmp/"
    message:
        "using fastq-dump to download SRA data files to {output.i1}"
    shell:"""
        parallel-fastq-dump -t 10  -O {params.DIR} --split-files --gzip -s {params.sample}
        mv {params.DIR}{params.sample}_1*.gz {params.DIR}{params.sample}_S1_L001_I1_001.fastq.gz
        mv {params.DIR}{params.sample}_2*.gz {params.DIR}{params.sample}_S1_L001_I2_001.fastq.gz
        mv {params.DIR}{params.sample}_3*.gz {params.DIR}{params.sample}_S1_L001_R1_001.fastq.gz
        mv {params.DIR}{params.sample}_4*.gz {params.DIR}{params.sample}_S1_L001_R2_001.fastq.gz
        """

rule count:
    input:
        i1 = WORKING_DIR + "fastq/{sample}_S1_L001_I1_001.fastq.gz",
        i2= WORKING_DIR + "fastq/{sample}_S1_L001_I2_001.fastq.gz",
        r1= WORKING_DIR + "fastq/{sample}_S1_L001_R1_001.fastq.gz",
        r2= WORKING_DIR + "fastq/{sample}_S1_L001_R2_001.fastq.gz"

    output: WORKING_DIR+"result_{sample}"
    threads: config["threads"]
    params:
       sample = "{sample}"
    shell:
        "cellranger count --id=result_{params.sample} \
                  --transcriptome=/home/mlei/ljscseq/medicago_truncatula/scrna-seq/Mtrun \
                  --fastqs=/home/mlei/ljscseq/medicago_truncatula/scrna-seq/fastq/ \
                  --sample={params.sample} "


