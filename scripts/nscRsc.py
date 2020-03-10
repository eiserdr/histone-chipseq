### python nscRsc.py <concatonatedPhantomFile> <outputFile>
### Written by Daina 
### Updated: 2020 Feb 13
### This python program will check files created by the cross-correlation program from phantompeakqualtools (run_spp.R). 
### The output file from run_spp.R calculates the NSC and RSC values for the -5 strand-shift. 
### However, we want the true strand shift around 160. This program will parse the phantomPeak file to recalculate the
### correct nsc and rsc files.  It will write to a file with the sample name, nsc, and rsc values.
### The phantompeak file has the following format:
### <Sample Name> <num Reads> <estimated fragment lengths> <correlation of est frag lengths> <phantompeak strandshift> <strandshift at lowest corr> <min corr> <NSC(col4/col8)> <RSC((col4-col8)/(col6-col8))> <QualityTag>
### MP2_27Ac_merged.tagAlign.gz	201770968	-5,170,350	0.682227116617166,0.668912440036005,0.658472037176196	70	0.654604	1500	0.649357	1.050619	6.264592	2
### More info at: https://github.com/kundajelab/phantompeakqualtools

###To use this program: python nscRsc.py <input_phantompeakfile> <new_output_filename>

import sys, getopt
def main(input, output):
	file = open(input, "r") #open phantomPeak file
	out = open(output,"a+") #open the output file
	out.write("Sample\tNSC\tRSC\n") #add a header to the output file
	list = file.readlines()
	for f in list:
		if f[0] != '#':							#ignores comments
			cols = f.split()					#split line into columns
			cc_peaks = cols[2].split(',')		#the peak values in the 3rd column are comma delimited
			cc_vals = cols[3].split(',')		#split the fourth column by comma
			if (int(cc_peaks[0]) <= 0):				#if the first peak is negative or 0, then it will recalculate the nsc and rsc values
				nsc = float(cc_vals[1]) / float(cols[7])
				rsc = (float(cc_vals[1]) - float(cols[7])) / (float(cols[5]) - float(cols[7]))
			else:
				nsc = float(cols[8])
				rsc = float(cols[9])
			print(cols[0])
			out.write(str(cols[0]) + '\t' + str(nsc) + '\t' + str(rsc) + '\n')  #output name, nsc, and rsc
	file.close()
	out.close()
		
if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])
