################################################################################
# TITLE:                    2020 RNAseq pipeline (Goodboy)
# AUTHOR:                   Ryan Rebernick
# DATE LAST MODIFIED:       04/19/2020
#
# FUNCTION:                 FastQC
#                           STAR alignment
#                           FeatureCounts (subread)
#                           MultiQC
#
# USES:
#                           - python 3.7.6
#                           - snakemake-minimal 5.4.0
#                           - fastQC v0.11.8
#                           - STAR STAR_2.6.1a_08-27
#                           - featureCounts v1.6.3
#                           - multiqc version 1.7
#
################################################################################

from os.path import join
import os
import glob


##############################################################
# Globals
##############################################################

# directories
LOG_DIR = config["OUT_DIR"] + "log/"
FQC_DIR = config["OUT_DIR"] + "fastqc/"
GENOME_DIR = config["GENOME_DIR"]
STAR_DIR = config["OUT_DIR"] + "star/"
FC_DIR = config["OUT_DIR"] + "featureCounts/"
MULTIQC_DIR = config["OUT_DIR"] + "multiqc/"

# samples/read names
READS = ["1", "2"]
SAMPLES = [os.path.basename(fname).split('.')[0] for fname in glob.glob(config["FQ_DIR"] + '*.R1.fastq.gz')]

# Genome FILES
FA_NAME = [os.path.basename(fname).split('.fa.gz')[0] for fname in glob.glob(GENOME_DIR + '/*.fa.gz')]
GTF_NAME = [os.path.basename(fname).split('.gtf.gz')[0] for fname in glob.glob(GENOME_DIR + '/*.gtf.gz')]
FA_GZ = expand(GENOME_DIR + "/{file}.fa.gz", file=FA_NAME)
FA = expand(GENOME_DIR + "/{file}.fa", file=FA_NAME)
GTF_GZ = expand(GENOME_DIR + "/{file}.gtf.gz", file=GTF_NAME)
GTF = expand(GENOME_DIR + "/{file}.gtf", file=GTF_NAME)

# FASTQ files
FQ = expand(config["FQ_DIR"] + "{sample}.R{read}.fastq.gz", sample=SAMPLES, read=READS)


##############################################################
# List of end point files needed
##############################################################

FQC = expand(FQC_DIR + "{sample}.R{read}_fastqc.html", sample=SAMPLES, read=READS)
GENOME = [GENOME_DIR + "/genomeParameters.txt"]
ALIGNED = expand(config["OUT_DIR"] + "star/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES)
FC = [FC_DIR + "counts.txt"]
MULTIQC = [MULTIQC_DIR + "multiqc_report.html"]


##############################################################
# Rules
##############################################################

# end files required
rule all:
    input: FQC + FA + GTF + GENOME + ALIGNED + FC + MULTIQC
    params:
    shell:
        """
        """


# FastQC
rule fastqc:
    input: FQ
    output: FQC
    params:
        fqc = FQC_DIR,
        log = LOG_DIR
    log: LOG_DIR + "FQC.log"
    shell: """ \
    mkdir -p {params.fqc}; \
    mkdir -p {params.log}; \
    fastqc -o {params.fqc} {input} \
    """


# UNZIP fa
rule unzip_fa:
    input: FA_GZ
    output: FA
    threads: 2
    params:
    shell:
        """
        gunzip -c {input} > {output};
        """


# UNZIP gtf
rule unzip_gtf:
    input: GTF_GZ
    output: GTF
    threads: 2
    params:
    shell:
        """
        gunzip -c {input} > {output};
        """


#Index the genome for STAR
rule starGenomeIndex:
    input:
        gd = GENOME_DIR,
        fa = FA,
        gtf = GTF
    output:
        index = GENOME
    params:
        readLength = config["READ_LENGTH"],
        gd = GENOME_DIR
    log: LOG_DIR + "genomeIndex.log"
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {input.gd} \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang {params.readLength};
        """

#STAR aligner
rule star:
    input:
        R1 = config["FQ_DIR"] + "{sample}.R1.fastq.gz",
        R2 = config["FQ_DIR"] + "{sample}.R2.fastq.gz",
        index = GENOME
    output:
        out = config["OUT_DIR"] + "star/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        ext = STAR_DIR + "{sample}_",
        gd = GENOME_DIR,
        star = STAR_DIR
    log: LOG_DIR + "star_{sample}.log"
    threads: 8
    shell:
        """
        mkdir -p {params.star}; \
        STAR \
        --genomeDir {params.gd} \
        --runThreadN {threads} \
        --readFilesIn {input.R1} {input.R2} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix {params.ext} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard
        """

rule featureCounts:
    input:
        bam = ALIGNED,
        gd = GENOME_DIR
    output: FC
    params:
        fc = FC_DIR,
        star = STAR_DIR
    log: LOG_DIR + "featureCounts.log"
    threads: 8
    shell:
        """
        mkdir -p {params.fc}; \
        featureCounts -p -t exon -g gene_id \
        -a {input.gd}/*.gtf \
        -o {params.fc}/counts.txt \
        {params.star}*Aligned.sortedByCoord.out.bam
        """


# MULTIQC
rule multiqc:
    input: FQC + GENOME + ALIGNED + FC
    output: MULTIQC
    params:
        outDir = MULTIQC_DIR,
        searchDir1 = config['OUT_DIR'],
        searchDir2 = config['FQ_DIR']
    threads: 8
    log: LOG_DIR + "multiqc.log"
    shell:
        """
        mkdir -p {params.outDir}; \
        multiqc -o {params.outDir} \
        {params.searchDir1} {params.searchDir2}
        """
