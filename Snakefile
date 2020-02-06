#sbatch --cpus-per-task=2 --time=2:00:00 snakemake.sh
#snakemake --unlock
#To remove everyting made by snakemake: snakemake all --delete-all-output
import sys

configfile: "config.yaml"
localrules: all

CONTROL = config["control"].split()		#there is only one control, i make it into a list so it can concatonate with the other lists
SAMPLES_BROAD = config["samples_broad"]
SAMPLES_NARROW = config["samples_narrow"]
SAMPLES_SICER = config["samples_sicer"]
ALL_CASES = SAMPLES_BROAD + SAMPLES_NARROW + SAMPLES_SICER
ALL_SAMPLES = CONTROL + ALL_CASES

ALL_BAM = expand("BamFiles/{sample}.sorted.bam", sample = ALL_SAMPLES)


ALL_PEAKS = expand("Callpeak/Broadpeak/{sample}_peaks.broadPeak", sample = SAMPLES_BROAD) + \
expand("Callpeak/Narrowpeak/{sample}_peaks.narrowPeak", sample = SAMPLES_NARROW) + \
expand("Callpeak/SICER/{sample}-W500-G1500-FDR0.01-island.bed", sample = SAMPLES_SICER)

ALL_SIGNAL = expand("Callpeak/Broadpeak/{sample}_{signal}.bw", sample = SAMPLES_BROAD, signal = ["ppois", "FE"]) + \
expand("Callpeak/Narrowpeak/{sample}_{signal}.bw", sample = SAMPLES_NARROW, signal = ["ppois", "FE"]) + \
expand("Callpeak/SICER/{sample}-W500-G1500-FDR0.01-islandfiltered-normalized.wig", sample = SAMPLES_SICER)

ALL_IDR = expand("IDR/BamFiles/{sample}_pr{rep}.bam", sample = ALL_CASES, rep =["1","2"])

rule all:
	input: ALL_BAM + ALL_PEAKS + ALL_SIGNAL + ALL_IDR

rule macs2_broadpeak:
	input:
		case="BamFiles/{sample}.sorted.bam",							#it knows that to only run it on samples specified by SAMPLES_BROAD
		ctrl = "BamFiles/" + config["control"] + ".sorted.bam"
	output: 
		"Callpeak/Broadpeak/{sample}_peaks.broadPeak", temp("Callpeak/Broadpeak/{sample}_treat_pileup.bdg"), temp("Callpeak/Broadpeak/{sample}_control_lambda.bdg"), temp("Callpeak/Broadpeak/{sample}_peaks.xls"), temp("Callpeak/Broadpeak/{sample}_peaks.gappedPeak")
	shell:
		"macs2 callpeak -t {input.case} -f AUTO -c {input.ctrl} -g hs -n {wildcards.sample} --outdir Callpeak/Broadpeak --broad --bw 150 --mfold 10 30 --bdg --nomodel --extsize 150  --SPMR"
#temp() will remove those files once the other rules are done using them. 

rule macs2_narrowpeak:
	input:
		case="BamFiles/{sample}.sorted.bam",
		ctrl = "BamFiles/" + config["control"] + ".sorted.bam"
	output: 
		"Callpeak/Narrowpeak/{sample}_peaks.narrowPeak", temp("Callpeak/Narrowpeak/{sample}_treat_pileup.bdg"), temp("Callpeak/Narrowpeak/{sample}_control_lambda.bdg"), temp("Callpeak/Narrowpeak/{sample}_peaks.xls"), temp("Callpeak/Narrowpeak/{sample}_summits.bed")
	shell:
		"macs2 callpeak -t {input.case} -f AUTO -c {input.ctrl} -g hs -n {wildcards.sample} --outdir Callpeak/Narrowpeak --bw 150 --mfold 10 30 --bdg --nomodel --extsize 150 --SPMR"
		#don't need {wildcards.sample}_peaks.xls, {wildcards.sample}_summits.bed
		
#sicer intakes, bed files. And the bed files need to be in separate directories, bc sicer's temporary directories will interfere with each other
rule bamtobed:
	input:
		"BamFiles/{sample}.sorted.bam"
	output:
		"BedFiles/{sample}.bed"
	shell:
		"bedtools bamtobed -i {input} > {output}"

#sicer2 won't run if I specify an out directory. So I have to make another rule to move the files to the right directory		
rule sicer:
	input:
		case="BedFiles/{sample}.bed",
		ctrl= "BedFiles/" + config["control"] + ".bed"
	output: 
		"{sample}-W500-G1500-FDR0.01-island.bed", "{sample}-W500-G1500-FDR0.01-islandfiltered-normalized.wig", temp("{sample}-W500-G1500-FDR0.01-islandfiltered.bed"), temp("{sample}-W500-G1500.scoreisland"), temp("{sample}-W500-normalized.wig"), temp("{sample}-W500-G1500-islands-summary")
	threads: 16
	shell:
		"sicer -t {input.case} -c {input.ctrl} -s hg19 -rt 1 -w 500 -f 150 -egf .8 -g 1500 -fdr 0.01 -cpu $(($SLURM_CPUS_PER_TASK/2)) --significant_reads"
rule mv_sicer:
	input:
		"{sample}-W500-G1500-FDR0.01-island.bed", "{sample}-W500-G1500-FDR0.01-islandfiltered-normalized.wig"
	output:
		"Callpeak/SICER/{sample}-W500-G1500-FDR0.01-island.bed", "Callpeak/SICER/{sample}-W500-G1500-FDR0.01-islandfiltered-normalized.wig"
	shell:
		"""
		mv {input} ./Callpeak/SICER
		"""
		#maybe remove the other files like "rm    
###Signal Files ####
#sicer already makes signal files, but I need to for macs2

rule bdgcmp:
	input:
		case="Callpeak/{dir}/{sample}_treat_pileup.bdg", #the wildcards.peak is to account for narrowpeak and broadpeak
		ctrl="Callpeak/{dir}/{sample}_control_lambda.bdg"
	output:
		temp("Callpeak/{dir}/{sample}_ppois.bdg"), temp("Callpeak/{dir}/{sample}_FE.bdg")
	shell:
		"macs2 bdgcmp -t {input.case} -c {input.ctrl} -m FE ppois --o-prefix Callpeak/{wildcards.dir}/{wildcards.sample}"

#module ucsc
rule bdgTobw:
	input: "Callpeak/{dir}/{sample}.bdg"
	output: "Callpeak/{dir}/{sample}.bw"
	threads: 12
	shell:
		"""
		sort -k1,1 -k2,2n {input} > {input}.sorted;
		bedGraphToBigWig {input}.sorted /fdb/genomebrowser/chrom.sizes/hg19/chrom.sizes {output}
		rm {input}.sorted
		"""
rule idr_pr:
	input: "BamFiles/{sample}.sorted.bam"
	output: "IDR/BamFiles/{sample}_pr1.bam", "IDR/BamFiles/{sample}_pr2.bam"
	shadow: "shallow"	#shallow will create a temporary directory. It will then delete files I don't need. It's nice here because I create a lot of intermediate files.
	shell:
		"""
		nlines=$(samtools view {input} | wc -l)
		nlines=$(((nlines+1)/2))
		samtools view -H {input} > header_{wildcards.sample}.sam
		samtools view {input} | shuf - | split -d -l $nlines - pr_{wildcards.sample}_
		cat header_{wildcards.sample}.sam pr_{wildcards.sample}_00 | samtools view -b - > IDR/BamFiles/{wildcards.sample}_pr1.bam
		cat header_{wildcards.sample}.sam pr_{wildcards.sample}_01 | samtools view -b - > IDR/BamFiles/{wildcards.sample}_pr2.bam
		"""
		
		
		
		
		
		
