## Script to convert UCSC/Ensembl gene annotation table dump into bed file
##
## How to generate input file
## In the UCSC genome browser, go to the table browser and select genome assembly
## (e.g., hg19). Select group "Genes and Gene Predictions", track "Ensembl
## Genes", and table ensGene. For output format, select "selected fields from
## primary and related tables". Enter an output file name, and click "get 
## output". From the primary table, select the following: name, chrom, strand, 
## txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, and
## name2. From "ensemblToGeneName", select name and value.
##

def run(inTab, outBed):
	import gzip
	import subprocess as sp

	if '.gz' in inTab:
		fin = gzip.open(inTab, 'rb')
	else:
		fin = open(inTab, 'r')
	ens = fin.read().split('\n')
	fin.close()

	header = ens[0]
	ens = ens[1:]
	if '' in ens:
		ens.remove('')

	fout = open(outBed+'.tmp', 'w')
	for line in ens:
		line = line.split('\t')

		## Process the exon data into proper format (relative to tss)
		txStartAbs = int(line[3]) # TSS
		exStartAbs = [int(n) for n in line[8].strip(',').split(',')] # Exon starts
		exStopAbs  = [int(n) for n in line[9].strip(',').split(',')] # Exon stops
		exSt = [ str(exStartAbs[i] - txStartAbs) for i in range(len(exStartAbs)) ]
		exLn = [ str(exStopAbs[i] - exStartAbs[i]) for i in range(len(exStartAbs)) ]

		## Adjust the order to fit bed12 standard + ENSG & Gene ID
		## the first '0' is for score, the second is for color (ucsc browser stuff)
		newline = '\t'.join([line[1], line[3], line[4], line[0], '0', line[2], line[5], 
							 line[6], '0', line[7], ','.join(exSt), ','.join(exLn), 
							 line[10], line[12]])
		fout.write(newline + '\n')
	fout.close()

	## Sort the bed file
	sp.call("less %s.tmp | sort -k1,1 -k2,2n | gzip > %s" % (outBed, outBed), 
			shell=True)
	sp.call("rm %s.tmp" % outBed, shell=True)
	

def main():
	import argparse

	parser = argparse.ArgumentParser()
	parser.add_argument("input_table", help="Input table to convert")
	parser.add_argument("output_bed", help="Name of output file")

	options = parser.parse_args()
	run(options.input_table, options.output_bed)

if __name__ == '__main__':
	main()
