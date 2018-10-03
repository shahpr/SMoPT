README file for Stochastic Model of Protein Translation (SMoPT)

*************************************************************************************

VERSION: SMoPT 2.0
DATE: June 28, 2015


ADDITIONAL INFO: see ./info folder for VERSION.txt, LICENSE.txt  and AUTHORS.txt

DESCRIPTION:
	The code simulates the dynamics of protein translation within an entire cell.
	The model keeps track of every ribosome, tRNA, and mRNA molecule.

	The model is described in: 
		Shah P, Ding Y, Niemczyk M, Kudla G, and Plotkin JB. 
		Rate-limiting steps in yeast protein translation. (2013) Cell 153 (7): 1589-1601

		Weinberg DE, Shah P, Eichhorn SW, Hussmann JA, Plotkin JB, and Bartel DP.
		Improved ribosome-footprint and mRNA measurements provide insights into dynamics and regulation of yeast translation. bioRxiv 10.1101/021501

	****When using this work please cite the above papers****

LICENSE: 
	Licensed under GNU General Public License (GPLv3). 
	See LICENSE.txt file for details.

BUILD:	
	Builds on Mac 10.8 and Ubuntu 12.04 machine using command: 
		gcc translation_v2.0.c -g -lm -lgsl -lgslcblas -mtune=generic -O3 -o SMoPT_v2

SYNOPSIS:

	./bin/SMoPT_V2 [options]

OPTIONS:

	-V <value>	Volume of the cell in m^3/s. The minimum volume of the cell is set to
			contain at least 1000 ribosomes and tRNAs.
			[DEFAULT]  -V 4.2E-17 (volume of yeast cell)

	-T[CHAR]		Specify times for various events

			-Tt:	Total simulation time in seconds.
				[DEFAULT]  -Tt 1500

			-Tb:	Burn-in/threshold time. Time spent by the cell to reach equilibrium.
				Only calculations after this time will be included in the analyses.
				[DEFAULT]  -Tb 1000

			-Th:	Time at which harringtonine is added to the cell.
				[DEFAULT]  -Th 1500

			-Tc:	Time at which cycloheximide is added to the cell.
				[DEFAULT]  -Tc 1500

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

	-J <FILE>	File containing the initial state of the system to begin simulations from.
			This file is an output of this simulation code containing '*_ribo_pos_*'
			[DEFAULT]  -C example/input/output_final_ribo_pos.out

	-x[INTEGER]	Specify parameters for cycloheximide action

			-x1 <value>	Probability of cycloheximide binding ribosomes
					[DEFAULT]  -x1 0

			-x2 <value>	Rate of cycloheximide dissociation from bound ribosomes
					[DEFAULT]  -x2 0 (irreversible binding)

			-y <value> 	Harringtonine rate for free ribosomes. This needs to be specified by the user.
					[DEFAULT]  -y 0

	-s <value>	Random number seed. *MUST SETUP*
			DEFAULT]  -s 1

	-O <prefix>	Specifies the prefix for the output files.
			[DEFAULT]  -O output


	-p[INTEGER]	Specify which output files to print

			-p1:	Generates a file of average elongation times
				of all codons.

			-p2:	Generates a file of total average elongation
				time of each gene.

			-p3:	Generates a file of average time between initiation
				events on mRNAs of each gene.

			-p4:	Generates a file of average number of free ribosomes,
				and free tRNAs of each type.

			-p5:	Generates a file of the final state of all mRNAs in a cell.
				It contains the positions of all bound ribosomes on mRNAs.

			-p6:	This generates two files:
					A file containing the amount of time wasted by stalled
					ribosomes on mRNAs of each gene.
																																											A file containing the time wasted by stalled ribosomes
					on each codon position of Gene 0 (first gene in the
					processed fasta file).
																																											This option significantly increases the total running time
					of the simulation. Use it with caution.

			-p7:	Generates a file of state of all mRNAs in a cell every second.
				This is similar to -p5 printed every second.
				This option significantly increases the total running time
				of the simulation. Use it with caution.

			-p8:	Generates a file of number of bound ribosomes at each position of a gene.

			-p9:	Generates a file of average (RPF and mRNA based) of bound ribosomes
				at each position of a gene.

			
BINARIES:
	*NIX and OSX:
	Precompiled binary for *NIX and OSX machines are included in bin/.  
	OSX binary was compiled with GCC 4.5.4 with the command: 
		gcc translation_v2.0.c -g -lm -lgsl -lgslcblas -mtune=generic -O3 -o SMoPT_v2
	Linux binary was compiled Ubuntu/Linaro with GCC 4.6.3 with the command: 
		gcc translation_v2.0.c -g -lm -lgsl -lgslcblas -mtune=generic -O3 -o SMoPT_v2_Linux

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
	./bin/SMoPT_v2 -Tt 1500 -Tb 1000 -R 200000 -t 3300000 -N 4839 -F example/input/S.cer.flash-freeze.genom -C example/input/S.cer.tRNA -s 1413 -O example/output/output -p1 -p2 -p3 -p4 -p5 -p6 -p7 -p8 -p9

UPDATES:
	Updates for this code can be found at the following website: 
	https://github.com/shahpr/SMoPT

BUGS:
	In case of any bugs or trouble with the code, send a mail to premal.shah@rutgers.edu
