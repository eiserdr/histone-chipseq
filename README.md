# Histone-ChIPSeq-Pipeline
## Updated 26 Feb 2020
This pipeline used in Amundadottir Lab for Pancreas cell line histone chip-seq data on the 2018(CoLo, HPDE, hTert, KP4, PATU, SU8686) and 2019(PANC-1 and MIAPaCa-2) cell lines. This is only for use on NIH's HPC biowulf

# Setup

Navigate to directory with bam files. (This pipeline won't align fastq files. You have to align them yourself).  
Then clone this repository into your directory:
    
    mkdir BamFiles
    mv *.bam BamFiles
    git clone https://github.com/eiserdr/histone-chipseq.git

    
Edit the config.yaml file with your sample names. They are organized by which peak caller you want to use.
**samples_broad** will use macs2 broadpeak caller. Used for H3K27Ac and H3K4me1.  
**samples_narrow** will use macs2 regular narrowpeak caller. Used for H3K4me3.  
**sicer** uses sicer, for very long peak regions. Used for H3K27me3.  

```
├ BamFiles                #This is an example of what input is expected. Please use this naming convention of *_1.sorted.bam
│ ├ MP2_In_1.sorted.bam
│ ├ MP2_In_2.sorted.bam
│ ├ MP2_27Ac_1.sorted.bam
│ ├ MP2_27Ac_2.sorted.bam
│ ├ MP2_27me3_1.sorted.bam
│ ├ MP2_27me3_2.sorted.bam
│ ├ MP2_4me3_1.sorted.bam
│ ├ MP2_4me3_2.sorted.bam
│ ├ MP2_4me1_1.sorted.bam
│ └ MP2_4me1_2.sorted.bam
│
├ histone-chipseq
│ ├ cluster.json
│ ├ config.yaml     # This is the file you need to edit for input. See below
│ ├ README.md
│ ├ rulegraph.png
│ ├ scripts
│ │ ├ nscRsc.py
│ │ ├ run_spp.R
│ │ └ wgEncodeDacMapabilityConsensusExcludable.bed  #These are blacklisted regions of hg19 from UCSC genomebrowser
│ ├ Snakefile
│ └ snakemake.sh    #This is the file to submit to the cluster
│
├ log               
│
└ out
  ├ BedFiles        #Bed files produced from bamfiles for SICER
  │ └ *.bed
  ├ Callpeak
  │ ├ Broadpeak
  │ │ ├ *_filt_peaks.broadPeak      #blacklist filtered peaks
  │ │ └ Signal
  │ │   ├ *_FE.bw       #Signal file based on fold enrichment
  │ │   └ *_ppois.bw    #Signal file based on poisson p-value
  │ ├ Narrowpeak
  │ │ ├ *_filt_peaks.narrowPeak
  │ │ └ Signal
  │ │   ├ *_FE.bw
  │ │   └ *_ppois.bw
  │ └ SICER
  │   ├ *_filt_W500-G1500-FDR0.01-island.bed
  │   └ Signal
  │     └ *-W500-G1500-FDR0.01-islandfiltered-normalized.wig
  ├ IDR
  │ ├ BamFiles
  │ │ ├ *_pr1.bam       
  │ │ └ *_pr2.bam
  │ ├ Callpeak
  │ │ └ *_p0.01_filt_peaks.sorted.broadPeak
  │ ├ IDR
  │ │ └ *_IDR0.05_filt_peaks.sorted.broadPeak
  │ └ peakcount.txt
  └ PhantomPeaks
    ├ *_xcorr.pdf
    ├ *_xcorr.txt
    └ *_NSC_RSC.txt

```
