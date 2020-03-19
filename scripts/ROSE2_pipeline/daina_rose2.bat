#!/bin/bash
#script to run ROSE2 program and produce Super enhancer bed files. To submit:
#sbatch --mem=10g daina_rose2.bat

module load samtools R bamliquidator python/2.7
###Need to index input and 27Ac bam files
samtools index ../../../BamFiles/Acinar_In_merged.sorted.bam
samtools index ../../../BamFiles/Acinar_27Ac_1.sorted.bam
samtools index ../../../BamFiles/Acinar_27Ac_2.sorted.bam
samtools index ../../../BamFiles/Acinar_27Ac_merged.sorted.bam

###make bed file. It just needs a .bed extension
mv ../../../out/Callpeak/Broadpeak/Acinar_27Ac_1_filt_peaks.broadPeak ../../../out/Callpeak/Broadpeak/Acinar_27Ac_1_filt_peaks.bed
mv ../../../out/Callpeak/Broadpeak/Acinar_27Ac_2_filt_peaks.broadPeak ../../../out/Callpeak/Broadpeak/Acinar_27Ac_2_filt_peaks.bed
mv ../../../out/Callpeak/Broadpeak/Acinar_27Ac_merged_filt_peaks.broadPeak ../../../out/Callpeak/Broadpeak/Acinar_27Ac_merged_filt_peaks.bed
###Run the main ROSE2 program
python ROSE2_main.py -g HG19 -i ../../../out/Callpeak/Broadpeak/Acinar_27Ac_1_filt_peaks.bed -r ../../../BamFiles/Acinar_27Ac_1.sorted.bam -c ../../../BamFiles/Acinar_In_merged.sorted.bam -o ../../../out/ROSE2 -t 2500 -s 12500
### Remove extraneous files
rm -r ../../../out/ROSE2/gff \
../../../out/ROSE2/mappedGFF
### Create a bed file of just super enhancers
num_lines=$(wc -l < ../../../out/ROSE2/Acinar_27Ac_1_filt_peaks_Enhancers_withSuper.bed)
super_lines=$(awk '/description=\"Super/ {print FNR}' ../../../out/ROSE2/Acinar_27Ac_1_filt_peaks_Enhancers_withSuper.bed)
cut_off=$(($num_lines - $super_lines))
tail -$cut_off ../../../out/ROSE2/Acinar_27Ac_1_filt_peaks_Enhancers_withSuper.bed > ../../../out/ROSE2/Acinar_27Ac_1_SuperEnhancers.bed

python ROSE2_main.py -g HG19 -i ../../../out/Callpeak/Broadpeak/Acinar_27Ac_2_filt_peaks.bed -r ../../../BamFiles/Acinar_27Ac_2.sorted.bam -c ../../../BamFiles/Acinar_In_merged.sorted.bam -o ../../../out/ROSE2 -t 2500 -s 12500
###Get rid of extraneous files
rm -r ../../../out/ROSE2/gff \
../../../out/ROSE2/mappedGFF
###Make a bed file of just super enhancers.
num_lines=$(wc -l < ../../../out/ROSE2/Acinar_27Ac_2_filt_peaks_Enhancers_withSuper.bed)
super_lines=$(awk '/description=\"Super/ {print FNR}' ../../../out/ROSE2/Acinar_27Ac_2_filt_peaks_Enhancers_withSuper.bed)
cut_off=$(($num_lines - $super_lines))
tail -$cut_off ../../../out/ROSE2/Acinar_27Ac_2_filt_peaks_Enhancers_withSuper.bed > ../../../out/ROSE2/Acinar_27Ac_2_SuperEnhancers.bed

python ROSE2_main.py -g HG19 -i ../../../out/Callpeak/Broadpeak/Acinar_27Ac_merged_filt_peaks.bed -r ../../../BamFiles/Acinar_27Ac_merged.sorted.bam -c ../../../BamFiles/Acinar_In_merged.sorted.bam -o ../../../out/ROSE2 -t 2500 -s 12500
rm -r ../../../out/ROSE2/gff \
../../../out/ROSE2/mappedGFF 
num_lines=$(wc -l < ../../../out/ROSE2/Acinar_27Ac_merged_filt_peaks_Enhancers_withSuper.bed)
super_lines=$(awk '/description=\"Super/ {print FNR}' ../../../out/ROSE2/Acinar_27Ac_merged_filt_peaks_Enhancers_withSuper.bed)
cut_off=$(($num_lines - $super_lines))
tail -$cut_off ../../../out/ROSE2/Acinar_27Ac_merged_filt_peaks_Enhancers_withSuper.bed > ../../../out/ROSE2/Acinar_27Ac_merged_SuperEnhancers.bed

###Remove a few move extra files
rm ../../../out/ROSE2/*12KB_STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt \
../../../out/ROSE2/*AllEnhancers* \
../../../out/ROSE2/*Enhancers_withSuperStretch.bed \
../../../out/ROSE2/*_SuperStretchEnhancers.table.txt \
../../../out/ROSE2/*Enhancers_withStretch.bed \
../../../out/ROSE2/*StretchEnhancers.table.txt \
../../../out/ROSE2/*Plot_points_stretch.png
###Change the extension back to broadPeak
mv ../../../out/Callpeak/Broadpeak/Acinar_27Ac_1_filt_peaks.bed ../../../out/Callpeak/Broadpeak/Acinar_27Ac_1_filt_peaks.broadPeak
mv ../../../out/Callpeak/Broadpeak/Acinar_27Ac_2_filt_peaks.bed ../../../out/Callpeak/Broadpeak/Acinar_27Ac_2_filt_peaks.broadPeak
mv ../../../out/Callpeak/Broadpeak/Acinar_27Ac_merged_filt_peaks.bed ../../../out/Callpeak/Broadpeak/Acinar_27Ac_merged_filt_peaks.broadPeak
