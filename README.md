# Histone-ChIPSeq-Pipeline
## Updated 26 Feb 2020
This pipeline used in Amundadottir Lab for Pancreas cell line histone chip-seq data on the 2018(CoLo, HPDE, hTert, KP4, PATU, SU8686) and 2019(PANC-1 and MIAPaCa-2) cell lines. This is only for use on NIH's HPC biowulf

# Setup

Navigate to directory with bam files. (This pipeline won't align fastq files. You have to align them yourself).  
Then clone this repository into your directory:
    
    mkdir BamFiles
    mv *.bam BamFiles
    git clone https://github.com/eiserdr/histone-chipseq.git
    mv BamFiles/ histone-chipseq/

    
Edit the config.yaml file with your sample names. They are organized by which peak caller you want to use.  Don't worry about replicates.
**samples_broad** will use macs2 broadpeak caller. Used for H3K27Ac and H3K4me1.  
**samples_narrow** will use macs2 regular narrowpeak caller. Used for H3K4me3.  
**sicer** uses sicer, for very long peak regions. Used for H3K27me3.  

```
working_directory
│ 
├ BamFiles
│ ├ MP2_In_T1.sorted.bam
│ ├ MP2_In_T2.sorted.bam
│ ├ MP2_In_merged.sorted.bam
│ ├ MP2_27Ac_T1.sorted.bam
│ ├ MP2_27Ac_T2.sorted.bam
│ ├ MP2_27Ac_merged.sorted.bam
│ ├ MP2_27me3_T1.sorted.bam
│ ├ MP2_27me3_T2.sorted.bam
│ ├ MP2_27me3_merged.sorted.bam
│ ├ MP2_4me3_T1.sorted.bam
│ ├ MP2_4me3_T2.sorted.bam
│ ├ MP2_4me3_merged.sorted.bam
│ ├ MP2_4me1_T1.sorted.bam
│ ├ MP2_4me1_T2.sorted.bam
│ └ MP2_4me1_merged.sorted.bam
│
├ cluster.json
│
├ config.yaml

├ README.md

├ rulegraph.png

├ scripts

├ Snakefile

├ snakemake.sh



```
