#!/bin/bash
#sbatch --mem=10g daina_rose2.bat

#module load samtools R bamliquidator python/2.7

#python ROSE2_main.py -g HG19 -i ../../Callpeak/Broadpeak/MP2_27Ac_merged_filt_peaks.bed -r ../../BamFiles/MP2_27Ac_merged.sorted.bam -c ../../BamFiles/MP2_In_merged.sorted.bam -o MP2_ROSE2 -t 2500 -s 12500
#moved some code out
#pipeline_template.py, pipeline_dfci.py, pipeline.py
#
#python ROSE2_main.py -g HG19 -i ../../Callpeak/Broadpeak/MP2_27Ac_merged_filt_peaks.bed -r ../../BamFiles/MP2_27Ac_merged.sorted.bam -c ../../BamFiles/MP2_In_merged.sorted.bam -o MP2_ROSE2_nopipeline -t 2500 -s 12500
#
#GPL16043.* dynamicEnhancer* commandline_template.py callBowtie* normalizeRNASeq.R python* RNA_SEQ_PIPELINE_README.txt tophatTemplate.sh extractGuides.py bamPlot*
#python ROSE2_main.py -g HG19 -i ../../Callpeak/Broadpeak/MP2_27Ac_merged_filt_peaks.bed -r ../../BamFiles/MP2_27Ac_merged.sorted.bam -c ../../BamFiles/MP2_In_merged.sorted.bam -o MP2_ROSE2_nofiles3 -t 2500 -s 12500
#

#rm bamToGFFExample.sh bamTableUpdate.py bamToGFF_turbo.py
#rm heatMapOrdered.R processGeckoBam.py makeBamMeta.py clusterEnhancer.*

###make bed file
#cp ../../../out/Callpeak/Broadpeak/KP4_27Ac_merged_filt_peaks.broadPeak ../../../out/Callpeak/Broadpeak/KP4_27Ac_merged_filt_peaks.bed
###index input and 27Ac bam files
#samtools index ../../../BamFiles/KP4_27Ac_merged.sorted.bam 
#samtools index ../../../BamFiles/KP4_In_merged.sorted.bam
#python ROSE2_main.py -g HG19 -i ../../../out/Callpeak/Broadpeak/KP4_27Ac_merged_filt_peaks.bed -r ../../../BamFiles/KP4_27Ac_merged.sorted.bam -c ../../../BamFiles/KP4_In_merged.sorted.bam -o ../../../out/ROSE2 -t 2500 -s 12500
rm -r ../../../out/ROSE2/gff \
../../../out/ROSE2/mappedGFF \
../../../out/ROSE2/*12KB_STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt \
../../../out/ROSE2/*AllEnhancers* \
../../../out/ROSE2/*Enhancers_withSuperStretch.bed \
../../../out/ROSE2/*_SuperStretchEnhancers.table.txt \
../../../out/ROSE2/*Enhancers_withStretch.bed \
../../../out/ROSE2/*StretchEnhancers.table.txt \
../../../out/ROSE2/*Plot_points_stretch.png