# Histone-ChIPSeq-Pipeline
## Updated 16 Mar 2020
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

If your samples are named KP4, here's a quick way to alter the config file.
```
sed -i 's/Sample/KP4/g' histone-chipseq/config.yaml 
```
The Bam Files should follow the naming convention shown below:

```
├ BamFiles                #This is an example of what input is expected. Please use this naming convention of *_1.sorted.bam
│ ├ Sample_In_1.sorted.bam
│ ├ Sample_In_2.sorted.bam
│ ├ Sample_27Ac_1.sorted.bam
│ ├ Sample_27Ac_2.sorted.bam
│ ├ Sample_27me3_1.sorted.bam
│ ├ Sample_27me3_2.sorted.bam
│ ├ Sample_4me3_1.sorted.bam
│ ├ Sample_4me3_2.sorted.bam
│ ├ Sample_4me1_1.sorted.bam
│ └ Sample_4me1_2.sorted.bam
│
├ histone-chipseq
│ ├ cluster.json
│ ├ config.yaml     # This is the file you need to edit for input. See below
│ ├ README.md
│ ├ rulegraph.png
│ ├ run_rose2.sh    #can be used to submit the daina_rose2.bat script
│ ├ scripts
│ │ ├ hg19-blacklist.v2.bed
│ │ ├ hg19.chrom.sizes
│ │ ├ hg38-blacklist.v2.bed
│ │ ├ hg38.chrom.sizes
│ │ ├ nscRsc.py
│ │ ├ run_spp.R
│ │ └ ROSE2_pipeline
│ │   ├ daina_rose2.bat     #script to run ROSE2
│ │   ├ ...
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
Callpeaks broadpeak is run with these parameters on macs2.2.5
```
macs2 callpeak -t {input.case} -f AUTO -c {input.ctrl} -g hs -n {wildcards.sample} --outdir out/Callpeak/Broadpeak --broad --bw 150 --mfold 10 30 --bdg --nomodel --extsize 150  --SPMR
```
Callpeaks narrowpeak is run with the same parameter, but without the --broad parameter
Signal Files are produced with macs2 bdgcmp:
```
macs2 bdgcmp -t {input.case} -c {input.ctrl} -m FE ppois
```

SICER is run on sicer2-1.0.0
```
sicer -t {input.case} -c {input.ctrl} -s hg19 -rt 1 -w 500 -f 150 -egf .8 -g 1500 -fdr 0.01 -cpu $(($SLURM_CPUS_PER_TASK/2)) --significant_reads
```
Phantom Peaks are produced by run_spp.R from phantompeakqualtools
The NSC RSC values produced by that is not using the right input, so I recalculate them useing nscRsc.py
```
Rscript histone-chipseq/scripts/run_spp.R  -c=*tagalign.gz -out={xcorr.text} -p=$SLURM_CPUS_PER_TASK -savp={xcorr.plot} -rf
python histone-chipseq/scripts/nscRsc.py {input} {output}
```

## ROSE2 
ROSE2 must be run separately. First, change the name of the files in the script daina_rose2.bat.  
Then you can submit run_rose.sh. 
```
#change the files names
sed -i 's/Sample/KP4/g' histone-chipseq/scripts/ROSE2_pipeline/daina_rose2.bat
#must submit the script in the same directory
cd histone-chipseq/scripts/ROSE2_pipeline
sbatch --mem=10g daina_rose2.bat
```
