################################################################################
# TITLE:                    2020 RNAseq pipeline (Goodboy)
# AUTHOR:                   Ryan Rebernick
# DATE LAST MODIFIED:       04/03/2020
#
# FUNCTION:                 FastQC
#                           STAR alignment
#                           FeatureCounts (subread)
#                           MultiQC
#
# USES:
#                           - python3.7-anaconda/2019.07
#                           - config/python_env/smk_env1.yaml
#                           - fastQC
#                           - STAR
#                           - subread
#                           - multiqc
#                           - config/python_env/multiqc_env1.yaml
#
# RUN: snakemake -s sundayMorning.smk --configfile sm_config.yaml --cores 1
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
FA_NAME = [os.path.basename(fname).split('.fa.gz')[0] for fname in glob.glob(GENOME_DIR + '*.fa.gz')]
GTF_NAME = [os.path.basename(fname).split('.gtf.gz')[0] for fname in glob.glob(GENOME_DIR + '*.gtf.gz')]

# FASTQ files
FQ = expand(config["FQ_DIR"] + "{sample}.R{read}.fastq.gz", sample=SAMPLES, read=READS)


##############################################################
# List of end point files needed
##############################################################

FQC = expand(FQC_DIR + "{sample}.R{read}_fastqc.html", sample=SAMPLES, read=READS)
GENOME = [GENOME_DIR + "genomeParameters.txt"]
ALIGNED = expand(config["OUT_DIR"] + "star/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES)
FC = [FC_DIR + "counts.txt"]
MULTIQC = [MULTIQC_DIR + "multiqc_report.html"]


##############################################################
# Rules
##############################################################

# end files required
rule all:
    input: FQC + GENOME + ALIGNED + FC
    params: a = 'start', b = SAMPLES, c = config["FQ_DIR"], d = 'stop'
    shell:
        """
        echo {params.a}
        echo {params.b}
        echo {params.c}
        echo {params.d}
        """


# FastQC
rule fastqc:
    input: FQ
    output: FQC
    params:
        fqc = FQC_DIR,
        log = LOG_DIR
    shell: """ \
    mkdir -p {params.fqc}; \
    mkdir -p {params.log}; \
    fastqc -o {params.fqc} {input} \
    """


# UNZIP fa
rule unzip_fa:
    input: expand(GENOME_DIR + "{file}.fa.gz", file=FA_NAME)
    output: expand(GENOME_DIR + "{file}.fa", file=FA_NAME)
    log: expand(config["OUT_DIR"] + "log/unzip.{file}.log", file=FA_NAME)
    threads: 2
    params:
    shell:
        """
        gunzip -c {input} > {output};
        """


# UNZIP gtf
rule unzip_gtf:
    input: expand(GENOME_DIR + "{file}.gtf.gz", file=GTF_NAME)
    output: expand(GENOME_DIR + "{file}.gtf", file=GTF_NAME)
    log: expand(LOG_DIR + "unzip.{file}.log", file=GTF_NAME)
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
        fa = rules.unzip_fa.output,
        gtf = rules.unzip_gtf.output
    output:
        index = GENOME
    params:
        readLength=config["READ_LENGTH"],
        gd = GENOME_DIR
    log: LOG_DIR + "genomeIndex.log"
    threads: 2
    shell:
        """
        mkdir -p {}
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {input.gd} \
        --genomeFastaFiles {input.fa}  \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang {params.readLength};
        """

#STAR aligner
rule star:
    input:
        files = FQ,
        index = GENOME
    output: ALIGNED
    params:
        ext = expand(STAR_DIR + "{sample}_", sample=SAMPLES),
        gd = GENOME_DIR,
        star = STAR_DIR
    threads: 6
    shell:
        """
        STAR \
        --genomeDir {params.gd} \
        --runThreadN {threads} \
        --readFilesIn {input.files} \
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
    threads: 2
    shell:
        """
        mkdir -p {params.fc}; \
        featureCounts -p -t exon -g gene_id \
        -a {input.gd}*.gtf \
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
    threads: 1
    conda: CONFIG_DIR + "multiqc_env1.yaml"
    shell:
        """
        mkdir -p {params.outDir}; \
        multiqc -o {params.outDir} \
        {params.searchDir1} {params.searchDir2}
        """
