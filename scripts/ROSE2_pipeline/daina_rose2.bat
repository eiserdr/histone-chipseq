#!/bin/bash
#sbatch --mem=10g daina_rose2.bat
###index input and 27Ac bam files
module load samtools R bamliquidator python/2.7
#samtools index ../../../BamFiles/KP4_In_merged.sorted.bam
#samtools index ../../../BamFiles/KP4_27Ac_1.sorted.bam
#samtools index ../../../BamFiles/KP4_27Ac_2.sorted.bam
#samtools index ../../../BamFiles/KP4_27Ac_merged.sorted.bam

###make bed file. It just needs a .bed extension
mv ../../../out/Callpeak/Broadpeak/KP4_27Ac_1_filt_peaks.broadPeak ../../../out/Callpeak/Broadpeak/KP4_27Ac_1_filt_peaks.bed
#mv ../../../out/Callpeak/Broadpeak/KP4_27Ac_2_filt_peaks.broadPeak ../../../out/Callpeak/Broadpeak/KP4_27Ac_2_filt_peaks.bed
#mv ../../../out/Callpeak/Broadpeak/KP4_27Ac_merged_filt_peaks.broadPeak ../../../out/Callpeak/Broadpeak/KP4_27Ac_merged_filt_peaks.bed


python ROSE2_main.py -g HG19 -i ../../../out/Callpeak/Broadpeak/KP4_27Ac_1_filt_peaks.bed -r ../../../BamFiles/KP4_27Ac_1.sorted.bam -c ../../../BamFiles/KP4_In_merged.sorted.bam -o ../../../out/ROSE2 -t 2500 -s 12500
rm -r ../../../out/ROSE2/gff \
../../../out/ROSE2/mappedGFF

#python ROSE2_main.py -g HG19 -i ../../../out/Callpeak/Broadpeak/KP4_27Ac_2_filt_peaks.bed -r ../../../BamFiles/KP4_27Ac_2.sorted.bam -c ../../../BamFiles/KP4_In_merged.sorted.bam -o ../../../out/ROSE2 -t 2500 -s 12500
#rm -r ../../../out/ROSE2/gff \
#../../../out/ROSE2/mappedGFF
#python ROSE2_main.py -g HG19 -i ../../../out/Callpeak/Broadpeak/KP4_27Ac_merged_filt_peaks.bed -r ../../../BamFiles/KP4_27Ac_merged.sorted.bam -c ../../../BamFiles/KP4_In_merged.sorted.bam -o ../../../out/ROSE2 -t 2500 -s 12500
#rm -r ../../../out/ROSE2/gff \
#../../../out/ROSE2/mappedGFF 
rm ../../../out/ROSE2/*12KB_STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt \
../../../out/ROSE2/*AllEnhancers* \
../../../out/ROSE2/*Enhancers_withSuperStretch.bed \
../../../out/ROSE2/*_SuperStretchEnhancers.table.txt \
../../../out/ROSE2/*Enhancers_withStretch.bed \
../../../out/ROSE2/*StretchEnhancers.table.txt \
../../../out/ROSE2/*Plot_points_stretch.png
mv ../../../out/Callpeak/Broadpeak/KP4_27Ac_1_filt_peaks.bed ../../../out/Callpeak/Broadpeak/KP4_27Ac_1_filt_peaks.broadPeak
#mv ../../../out/Callpeak/Broadpeak/KP4_27Ac_2_filt_peaks.bed ../../../out/Callpeak/Broadpeak/KP4_27Ac_2_filt_peaks.broaPeak
#mv ../../../out/Callpeak/Broadpeak/KP4_27Ac_merged_filt_peaks.bed ../../../out/Callpeak/Broadpeak/KP4_27Ac_merged_filt_peaks.broadPeak