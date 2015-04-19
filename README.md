README file for Stochastic Model of Protein Translation (SMoPT)

*************************************************************************************

VERSION: SMoPT 1.0
DATE: May 06, 2013


ADDITIONAL INFO: see ./info folder for VERSION.txt, LICENSE.txt  and AUTHORS.txt

DESCRIPTION:
	The code simulates the dynamics of protein translation within an entire cell.
	The model keeps track of every ribosome, tRNA, and mRNA molecule.

	The model is described in: 
	Shah P, Ding Y, Niemczyk M, Kudla G, and Plotkin JB. 
	Rate-limiting steps in yeast protein translation. (2013) Cell 153 (7): 1589-1601

	****When using this work please cite the above paper****

LICENSE: 
	Licensed under GNU General Public License (GPLv3). 
	See LICENSE.txt file for details.

BUILD:	
	Builds on Mac 10.8 and Ubuntu 12.04 machine using command: 
		gcc translation.c -g -lm -mtune=generic -O3 -o SMoPT

SYNOPSIS:

	./bin/SMoPT [options]

OPTIONS:

	-V <value>	Volume of the cell in m^3/s. The minimum volume of the cell is set to 
			contain at least 1000 ribosomes and tRNAs.
			[DEFAULT]  -V 4.2E-17 (volume of yeast cell)

	-T <value>	Total simulation time in seconds.
			[DEFAULT]  -T 1500

	-H <value>	Burn-in/threshold time. Time spent by the cell to reach equilibrium.
			Only calculations after this time will be included in the analyses.
			[DEFAULT]  -H 1000

	-R <value>	Total number of ribosomes in the cell.
			[DEFAULT]  -R 200000

	-t <value>	Total number of tRNAs in the cell.
			[DEFAULT]  -t 3300000

	-N <value>	Total number of genes. This needs to be specified by the user.
			[DEFAULT]  -N 1

	-F <FILE>	File containing processed fasta file into a numeric sequence.
			This file is an output of the code utilities/convert.fasta.to.genom.pl
			It contains the information regarding initiation probability, mRNA
			abundance and codon sequence of each gene.
			[DEFAULT]  -F example/input/S.cer.genom

	-C <FILE>	File containing the information about codon, tRNA, tRNA abundance and wobble.
			This file is an output of the code utilities/create.Scer.cod.anticod.numeric.pl
			[DEFAULT]  -C example/input/S.cer.tRNA

	-s <value>	Random number seed.
			[DEFAULT]  -s 1

	-O <prefix>	Specifies the prefix for the output files.
			[DEFAULT] -O output


	-p[INTEGER]	Specify which output files to print
				
			-p1: Generates a file of average elongation times 
			     of all codons.

			-p2: Generates a file of total average elongation
			     time of each gene.

			-p3: Generates a file of average time between initiation
			     events on mRNAs of each gene.

			-p4: Generates a file of average number of free ribosomes,
			     and free tRNAs of each type.

			-p5: Generates a file of the final state of all mRNAs in a cell.
			     It contains the poistions of all bound ribosomes on mRNAs.

			-p6: This generates two files:
			     A file containing the amount of time wasted by stalled 
			     ribosomes on mRNAs of each gene.
			     A file containing the time wasted by stalled ribosomes 
			     on each codon position of Gene 0 (first gene in the 
			     processed fasta file).
			     This option significantly increases the total running time
			     of the simulation. Use it with caution.

			-p7: Generates a file of state of all mRNAs in a cell every second.
			     This is similar to -p5 printed every second.
			     This option significantly increases the total running time
			     of the simulation. Use it with caution.

			
BINARIES:
	*NIX and OSX:
	Precompiled binary for *NIX machines is included as bin/SMoPT.  
	This binary was compiled on a linux machine GCC 4.2 with the command: 
		gcc translation.c -g -lm -mtune=generic -O3 -o SMoPT
	The binary ideally should run on all i386 and x86_64 machines running linux or OSX.
	If problems are encountered we encourage you to try recompiling the code for your local 
	machine before contacting the authors.


COMPILERS:
	The code can be recompiled from source and optimized for the hardware of the 
	machine it will be run on.
	
	The source code has been successfully compiled with the following compilers.
	Mac:		gcc
	Linux:		gcc


EXAMPLES:
	(Fast version)	./bin/SMoPT -T 1500 -H 1000 -R 200000 -t 3300000 -N 3795 -F example/input/S.cer.genom -C example/input/S.cer.tRNA -s 1413 -O example/output/output -p1 -p2 -p3 -p4 -p5
	(Slow version)	./bin/SMoPT -T 1500 -H 1000 -R 200000 -t 3300000 -N 3795 -F example/input/S.cer.genom -C example/input/S.cer.tRNA -s 1413 -O example/output/output -p1 -p2 -p3 -p4 -p5 -p6 -p7

UPDATES:
	Updates for this code can be found at the following website: 
	http://mathbio.sas.upenn.edu/shah-cell-2013-code.tar.gz

BUGS:
	In case of any bugs or trouble with the code, send a mail to shahpr@sas.upenn.edu
