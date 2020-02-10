#sbatch --cpus-per-task=2 --time=2:00:00 snakemake.sh
#snakemake --unlock
#To remove everyting made by snakemake: snakemake all --delete-all-output
import sys

configfile: "config.yaml"
localrules: all

CONTROL = config["control"]		#there is only one control, note that the CONTROL variable is not a list.
SAMPLES_BROAD = config["samples_broad"]
SAMPLES_NARROW = config["samples_narrow"]
SAMPLES_SICER = config["samples_sicer"]
ALL_CASES = SAMPLES_BROAD + SAMPLES_NARROW + SAMPLES_SICER
ALL_SAMPLES = CONTROL.split() + ALL_CASES

ALL_BAM = expand("BamFiles/{sample}_{rep}.sorted.bam", sample = ALL_SAMPLES, rep = ["1","2","merged"])

ALL_PEAKS = expand("Callpeak/Broadpeak/{sample}_{rep}_filt_peaks.broadPeak", sample = SAMPLES_BROAD, rep = ["1","2","merged"]) + \
expand("Callpeak/Narrowpeak/{sample}_{rep}_filt_peaks.narrowPeak", sample = SAMPLES_NARROW, rep = ["1","2","merged"]) + \
expand("Callpeak/SICER/{sample}_{rep}-W500-G1500-FDR0.01-island.bed", sample = SAMPLES_SICER, rep = ["1","2","merged"])

ALL_SIGNAL = expand("Callpeak/Broadpeak/{sample}_{rep}_{signal}.bw", sample = SAMPLES_BROAD, rep = ["1","2","merged"], signal = ["ppois", "FE"]) + \
expand("Callpeak/Narrowpeak/{sample}_{rep}_{signal}.bw", sample = SAMPLES_NARROW, rep = ["1","2","merged"], signal = ["ppois", "FE"]) + \
expand("Callpeak/SICER/{sample}_{rep}-W500-G1500-FDR0.01-islandfiltered-normalized.wig", sample = SAMPLES_SICER, rep = ["1","2","merged"])

ALL_IDR = expand("IDR/BamFiles/{sample}_{rep}_pr{pr_rep}.bam", sample = SAMPLES_BROAD + SAMPLES_NARROW, rep = ["1","2","merged"], pr_rep =["1","2"]) + \
expand("IDR/Callpeak/{sample}_{rep}_p0.01_filt_peaks.sorted.narrowPeak", sample = SAMPLES_NARROW, rep = ["1","2","merged"]) + \
expand("IDR/Callpeak/{sample}_{rep}_p0.01_filt_peaks.sorted.broadPeak", sample = SAMPLES_BROAD, rep = ["1","2","merged"]) + \
expand("IDR/Callpeak/{sample}_{rep}_pr{pr_rep}_p0.01_filt_peaks.sorted.broadPeak", sample = SAMPLES_BROAD, rep = ["1","2","merged"], pr_rep = ["1","2"]) + \
expand("IDR/Callpeak/{sample}_{rep}_pr{pr_rep}_p0.01_filt_peaks.sorted.narrowPeak", sample = SAMPLES_NARROW, rep = ["1","2","merged"], pr_rep = ["1","2"])

rule all:
	input: ALL_BAM + ALL_PEAKS + ALL_SIGNAL + ALL_IDR

rule merge:
	input: "BamFiles/{sample}_1.sorted.bam", "BamFiles/{sample}_2.sorted.bam"
	output: "BamFiles/{sample}_merged.sorted.bam"
	shell:
		"samtools merge {output} {input}"

rule macs2_broad:
	input:
		case="BamFiles/{sample}.sorted.bam",							#it knows that to only run it on samples specified by SAMPLES_BROAD
		ctrl = "BamFiles/" + CONTROL + "_merged.sorted.bam"			#I only use the merged input here
	output: 
		"Callpeak/Broadpeak/{sample}_peaks.broadPeak", 
		temp("Callpeak/Broadpeak/{sample}_treat_pileup.bdg"), 
		temp("Callpeak/Broadpeak/{sample}_control_lambda.bdg"), 
		temp("Callpeak/Broadpeak/{sample}_peaks.xls"), 
		temp("Callpeak/Broadpeak/{sample}_peaks.gappedPeak")
	shell:
		"macs2 callpeak -t {input.case} -f AUTO -c {input.ctrl} -g hs -n {wildcards.sample} --outdir Callpeak/Broadpeak --broad --bw 150 --mfold 10 30 --bdg --nomodel --extsize 150  --SPMR"
#temp() will remove those files once the other rules are done using them. 

rule macs2_narrow:
	input:
		case="BamFiles/{sample}.sorted.bam",
		ctrl = "BamFiles/" + CONTROL + "_merged.sorted.bam"
	output: 
		"Callpeak/Narrowpeak/{sample}_peaks.narrowPeak", 
		temp("Callpeak/Narrowpeak/{sample}_treat_pileup.bdg"), 
		temp("Callpeak/Narrowpeak/{sample}_control_lambda.bdg"), 
		temp("Callpeak/Narrowpeak/{sample}_peaks.xls"), 
		temp("Callpeak/Narrowpeak/{sample}_summits.bed")
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
		ctrl= "BedFiles/" + CONTROL + "_merged.bed" #Only use merged input file
	output: 
		temp("{sample}-W500-G1500-FDR0.01-island.bed"), #This peak file is temp because it will be blacklist filtered
		"{sample}-W500-G1500-FDR0.01-islandfiltered-normalized.wig", 
		temp("{sample}-W500-G1500-FDR0.01-islandfiltered.bed"), 
		temp("{sample}-W500-G1500.scoreisland"), 
		temp("{sample}-W500-normalized.wig"), 
		temp("{sample}-W500-G1500-islands-summary")
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
		case="Callpeak/{dir}/{sample}_treat_pileup.bdg", 
		ctrl="Callpeak/{dir}/{sample}_control_lambda.bdg"
	output:
		temp("Callpeak/{dir}/{sample}_ppois.bdg"), temp("Callpeak/{dir}/{sample}_FE.bdg")  #These are temporary files that will be made into bigwig files.
	shell:
		"macs2 bdgcmp -t {input.case} -c {input.ctrl} -m FE ppois --o-prefix Callpeak/{wildcards.dir}/{wildcards.sample}"

#module ucsc
#couldn't figure out a way to write the bw output to the Signal file.  So I had to make two mv_bw rules
rule bdgTobw:
	input: "{dir}/{sample}.bdg"
	output: "{dir}/{sample}.bw"
	threads: 12
	shell:
		"""
		sort -k1,1 -k2,2n {input} > {input}.sorted
		bedGraphToBigWig {input}.sorted /fdb/genomebrowser/chrom.sizes/hg19/chrom.sizes {output}
		rm {input}.sorted
		"""

### IDR ###
		
rule make_pseudorep_idr:
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

##call peaks at low threshold of p0.01		
rule macs2_broad_idr:
	input: 
		case="BamFiles/{sample}_{rep}.sorted.bam",
		ctrl="BamFiles/" + CONTROL + "_{rep}.sorted.bam"
	output: 
		temp("IDR/Callpeak/{sample}_{rep}_p0.01_peaks.broadPeak"), #temp bc I keep sorted and filtered peakfile
		temp("IDR/Callpeak/{sample}_{rep}_p0.01_peaks.xls"), 
		temp("IDR/Callpeak/{sample}_{rep}_p0.01_peaks.gappedPeak")
	shell:
		"macs2 callpeak -p 0.01 -t {input.case} -f AUTO -c {input.ctrl} -g hs -n {wildcards.sample}_{wildcards.rep}_p0.01 --outdir IDR/Callpeak --broad --bw 150 --extsize 150 --nomodel"

rule macs2_narrow_idr:
	input: 
		case="BamFiles/{sample}_{rep}.sorted.bam",
		ctrl="BamFiles/" + CONTROL + "_{rep}.sorted.bam"
	output: 
		temp("IDR/Callpeak/{sample}_{rep}_p0.01_peaks.narrowPeak"), 
		temp("IDR/Callpeak/{sample}_{rep}_p0.01_peaks.xls"), 
		temp("IDR/Callpeak/{sample}_{rep}_p0.01_summits.bed")
	shell:
		"macs2 callpeak -p 0.01 -t {input.case} -f AUTO -c {input.ctrl} -g hs -n {wildcards.sample}_{wildcards.rep}_p0.01 --outdir IDR/Callpeak --bw 150 --extsize 150 --nomodel"
		
rule macs2_broad_pr_idr:
	input: 
		case="IDR/BamFiles/{sample}_{rep}_pr{pr_rep}.bam",
		ctrl="BamFiles/" + CONTROL + "_{rep}.sorted.bam"
	output: 
		temp("IDR/Callpeak/{sample}_{rep}_pr{pr_rep}_p0.01_peaks.broadPeak"),  
		temp("IDR/Callpeak/{sample}_{rep}_pr{pr_rep}_p0.01_peaks.xls"), 
		temp("IDR/Callpeak/{sample}_{rep}_pr{pr_rep}_p0.01_peaks.gappedPeak")
	shell:
		"macs2 callpeak -p 0.01 -t {input.case} -f AUTO -c {input.ctrl} -g hs -n {wildcards.sample}_{wildcards.rep}_pr{wildcards.pr_rep}_p0.01 --outdir IDR/Callpeak --broad --bw 150 --extsize 150 --nomodel"		
rule macs2_narrow_pr_idr:
	input: 
		case="IDR/BamFiles/{sample}_{rep}_pr{pr_rep}.bam",
		ctrl="BamFiles/" + CONTROL + "_{rep}.sorted.bam"
	output: 
		temp("IDR/Callpeak/{sample}_{rep}_pr{pr_rep}_p0.01_peaks.narrowPeak"), 
		temp("IDR/Callpeak/{sample}_{rep}_pr{pr_rep}_p0.01_peaks.xls"), 
		temp("IDR/Callpeak/{sample}_{rep}_pr{pr_rep}_p0.01_summits.bed")
	shell:
		"macs2 callpeak -p 0.01 -t {input.case} -f AUTO -c {input.ctrl} -g hs -n {wildcards.sample}_{wildcards.rep}_pr{wildcards.pr_rep}_p0.01 --outdir IDR/Callpeak --bw 150 --extsize 150 --nomodel"	

#Blacklist and sort_peaks both run on the idr peaks.  Their order doesn't matter. But snakemake calls it an ambiguous rule exception because it doesn't know which one to run first. I just tell snakemake to run blacklist first
ruleorder: blacklist > sort_peaks

#For all peak files. This removes blacklisted regions in hg19.
rule blacklist:
	input:
		"{Dir}/{sample}_{peaks}"
	output:
		"{Dir}/{sample}_filt_{peaks}"
	shell:
		"""
		bedtools intersect -v -a {input} -b Files/wgEncodeDacMapabilityConsensusExcludable.bed > {output}
		rm {input}
		"""	

#have to sort peaks by significance after blacklisting. I specify "filt" so that the blacklist occurs before sorting.  I don't think the order matters too much, but snakemake needs it.
rule sort_peaks:
	input:
		"{dir}/{sample}.{peak}"
	output:
		"{dir}/{sample}.sorted.{peak}"
	shell:
		"""
		sort -k8,8nr {input} > {output}
		rm {input}
		"""
