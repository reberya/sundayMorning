# sundayMorning
Basic paired end RNA-seq Analysis Pipeline availible for use on any system with sufficient memory capabilities. 



## Getting Started

You will need access to a computing system with: 
  - at least 16GB of memory with more potentially required depending on sample size.
  - [Singularity](https://sylabs.io/docs/) installed



### Prerequisites
Begin by cloning sundayMorning into your directory of choice. This will contain the snakemake file, snakemake configuaration file, and recipe file used to generate the singularity image. 
```
git clone <LINK_TO_SUNDAYMORNING>
```



### Getting the image file 

Next, navigate to sundayMorning/code and download the image file from sylabs cloud
```
singularity pull --arch amd64 library://reberya/default/sundaymorning:v1.0.0
```


### File modification
Several files must be modified to suit your specific installation of sundayMorning. Additionally a reference genome for alignment and feature counting is required.


#### Reference genome
Make your genome directory and get required files. Only store one .gtf.gz file and one .fa.gz in each genome folder if you'll be using multiple annotations. Storing multiple files will cause the program to throw errors. 
```
mkdir <YOUR_GENOME_FOLDER>
cd <YOUR_GENOME_FOLDER>
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/GRCm38.primary_assembly.genome.fa.gz
```
#### Read length
Change the config/sm_config.yaml file's READ_LENGTH parameter to equal the read length of your fastq files. 
```
READ_LENGTH = <YOUR_READ_LENGTH>
```

If you do not know this offhand you can run fastqc on one of your files from the sundaymorning:v1.0.0 shell and then check the report it generates.

```
singularity shell <.../.../sundayMorning/code/sundaymorning_v1.0.0.sif>

fastqc <FASTQ_FILE>
```


### Data organization
In order for sundayMorning to run without error on the availible fastq files, the folders containing your data must be structured appropriately. Appropriate structure consists of a single folder 'FQ' which contains your .fastq.gz files in the following format: <SAMPLE>.R<#>.fastq.gz where SAMPLE is your sample name and # is the read number [1,2]. Example file names would be:
  
sample1.R1.fastq.gz
sample1.R2.fastq.gz

or...

14524_M.R1.fastq.gz
14524_M.R2.fastq.gz



## Running sundayMorning
Once you have created your image file, downloaded the gtf/fasta annotation files, and ensured your data is organized appropriately you are ready to start your analysis. Begin with a dryrun to ensure snakemake can build properly and there are no errors

```
singularity exec --bind <YOUR/FASTQ/DIRECTORY>:/in --bind <YOUR/CHOSEN/OUTPUT/DIRECTORY>:/out --bind <YOUR/CHOSEN/GENOME/DIRECTORY>:/genome <.../.../sundayMorning/code/sundaymorning_v1.0.0.sif> snakemake -s <.../.../sundayMorning/code/sundayMorning.smk> --configfile <.../.../sundayMorning/code/sm_config.yaml> --dryrun
```

Then, simply run the following command after changing the directories in <> to suit match your local filesystem and cores to match the number of cores availible to you. 
```
singularity exec --bind <YOUR/FASTQ/DIRECTORY>:/in --bind <YOUR/CHOSEN/OUTPUT/DIRECTORY>:/out --bind <YOUR/CHOSEN/GENOME/DIRECTORY>:/genome <.../.../sundayMorning/code/sundaymorning_v1.0.0.sif> snakemake -s <.../.../sundayMorning/code/sundayMorning.smk> --configfile <.../.../sundayMorning/code/sm_config.yaml> --cores 8
```

If you are working off of a cluster with slurm capabilities such as the [University of Michigan Greatlakes](https://arc-ts.umich.edu/greatlakes/) cluster, I first recommend allocating yourself the proper resources using:
```
srun --nodes 8 --account=<YOUR_ACCOUNT> --time 16:00:00 --ntasks-per-node=4 --mem-per-cpu=16GB --pty /bin/bash
```


### Output
sundayMorning will create an output directory in <YOUR/CHOSEN/OUTPUT/DIRECTORY> containing the output for each tool stored in a seperate folder. The output/featureCounts/counts.txt file may be used for gene expression with programs like [DESEQ2](https://github.com/mikelove/DESeq2) or [EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html).

I recommend first starting with the MultiQC output to check sample quality and success of the run. It is located in output/multiqc/multiqc_report.html



## Built With

* [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - FASTQC analysis
* [STAR](https://github.com/alexdobin/STAR) - Alignment of Fastqc files
* [SUBREAD](http://subread.sourceforge.net/) - Feature Counting
* [MULTIQC](https://multiqc.info/) - Overall QC (compiles samples)



## Acknowledgments

Shout out to [Bennet Fauber](https://arc-ts.umich.edu/staff-member/bennet-fauber/) for his help configuring singularity. The sm.def file is based upon his code. Also thanks to [Kelly Sovacool](https://github.com/kelly-sovacool) for her help setting up snakemake on greatlakes server initially. The initial configuration files and base snakemake file are based off of her mwe_hpc_snakemake repository. 


