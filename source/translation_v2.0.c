// Code to simulate the translation process based on explicit diffusion and mass action properties of tRNAs, ribosome and mRNAs.

/*
To compile and run the code run:

	gcc translation.c -g -lm -mtune=generic -O3 -o SMoPT
	./SMoPT

*/
// Declaring Header Files
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h> 
#include <string.h>
#include <gsl/gsl_rng.h>

// Fixed parameters
#define MAX_GENES 5000				// Maximum number of genes supported
#define MAX_GENE_LEN_ALLW 5000		// Maximum allowed gene length in codons
#define MAX_TIME 2400000			// Maximum time of the simulation
#define char_len_tRNA 1.5e-8		// Characteristic length of tRNA
#define char_len_ribo 3e-8			// Characteristic length of ribosome
#define char_time_tRNA 5.719e-4		// Characteristic time of movement for tRNA (4.45e-7 * 1285.1)
#define char_time_ribo 5e-4			// Characteristic time of movement for ribosome

// Default global variables
int seed = 0;						// Seed for RNG
int n_genes = 1;					// Number of genes
int tot_ribo = 2e5;					// Total ribosomes
int tot_tRNA = 3.3e6;				// Total tRNAs
int tot_mRNA = 0;					// Total mRNA (initialization)
double tot_space = 4.2e-17;			// Total space within a cell
double avail_space_t = 1.24e7;		// Available space for tRNAs
double avail_space_r = 1.56e6;		// Available space for ribosomes
double tot_time = 1500;				// Total time for simulation
double thresh_time = 1000;			// Threshold time for analysis of e_times
double chx_rb_prob = 0.0;			// Cycloheximide action probability for bound ribosomes
double chx_ds_rate = 0.0;			// Cycloheximide dissociation rate from bound ribosomes
double harr_rf_rate = 0.0;			// Harringtonine action rate for free ribosomes
double chx_time = 1500;				// Threshold time for action of cycloheximide
double harr_time = 1500;			// Threshold time for action of harringtonine

// Run options
int printOpt[9] = {0,0,0,0,0,0,0,0,0};
char out_prefix[150] = "output";	// Prefix for output file names
char out_file[150];
char fasta_file[150] = "example/input/S.cer.genom";
char code_file[150] = "example/input/S.cer.tRNA";
char state_file[150] = "";



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Create ribosome and mRNA based structures
typedef struct
{	int mRNA;					// Bound to which mRNA
	int pos;					// Position on mRNA
	double t_trans_ini;			// Time of translation initiation
	double t_elong_ini;			// Time of arrival at current codon
	int elng_cod_list;			// The id of the list of elongatable codons
	int elng_pos_list;			// Position in the list of elongatable codons
	int inhbtr_bound;			// Is the ribosome bound with CHX or Harr?
} ribosome;

typedef struct
{	int gene;					// Gene id
	int ini_n;					// Number of initiation events
	int trans_n;				// Number of translation events
	double last_ini;			// Time of last initiation event
	double avg_time_to_ini;		// Average time to initiation
	double avg_time_to_trans;	// Average translation time
} transcript;

typedef struct
{	int seq[MAX_GENE_LEN_ALLW];	// Codon sequence of the gene
	int len;					// Length of the gene
	int exp;					// Gene expression level
	double ini_prob;			// Initiation probability of the mRNA
} gene;

typedef struct
{	char codon[4];				// Codon
	int tid;					// tRNA id
	int gcn;					// tRNA gene copy number
	double wobble;				// Wobble parmeter
} trna;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Reading the processed sequence file
int Read_FASTA_File(char *filename, gene *Gene)
{	FILE *fh;
	int c1=0,c2=0;
	char curr_char;
	
	fh=fopen(filename, "r");

	if(!fh)					// Check if file exists
	{	printf("\nModified FASTA/Sequence File Doesn't Exist\n");
		fflush(stdout);
		Help_out();
		exit(1);
	}

	fscanf(fh,"%lf",&Gene[c1].ini_prob);

	do
	{	fscanf(fh,"%d",&Gene[c1].exp);
		
		c2 = 0;
		do
		{	fscanf(fh,"%d",&Gene[c1].seq[c2]);
			c2++;
			curr_char = fgetc(fh);
		}while(curr_char != '\n');
		Gene[c1].len = c2;

		c1++;
	}while(fscanf(fh,"%lf",&Gene[c1].ini_prob) ==1);
}



// Read in tRNA file
int Read_tRNA_File(char *filename, trna *cTRNA)
{	FILE *fh;
	int c1=0,c2=0;
	char curr_char;
	
	fh=fopen(filename, "r");
	
	if(!fh)					// Check if file exists
	{	printf("\ntRNA File Doesn't Exist\n");
		fflush(stdout);
		Help_out();
		exit(1);
	}
	c1 = 0;
	c2 = 0;
	while(fscanf(fh,"%s",&cTRNA[c1].codon)==1)
	{	fscanf(fh,"%d%d%lf",&cTRNA[c1].tid,&cTRNA[c1].gcn,&cTRNA[c1].wobble);
		c1++;
	}
}



// Reading the state of the system
int Read_STATE_File(char *filename, int **R_grid)
{	FILE *fh;
	int c1=0,c2=0;
	char curr_char;
	
	fh=fopen(filename, "r");
	
	if(!fh)					// Check if file exists
	{	printf("\nState File Doesn't Exist\n");
		fflush(stdout);
		Help_out();
		exit(1);
	}
	
	c1 = 0;
	c2 = -1;
	do
	{	c2++;
		do
		{	fscanf(fh,"%d",&R_grid[c1][c2]);
			c2++;
			curr_char = fgetc(fh);
		}while(curr_char != '\n');
		c1++;
		c2 = 0;
	}while(fscanf(fh,"%d",&R_grid[c1][c2]) ==1);
}


// Help output
int Help_out()
{	printf("\nUsage:\n");
	printf("\t./bin/SMoPT [options]\n\n");
	printf("Options:\n");
	printf("\n\t-V <value>	Volume of the cell in m^3/s. The minimum volume of the cell is set to\n");
	printf("\t\t\tcontain at least 1000 ribosomes and tRNAs.\n");
	printf("\t\t\t[DEFAULT]  -V 4.2E-17 (volume of yeast cell)\n");
	printf("\n");
	printf("\t-T[CHAR]	Specify times for various events\n");
	printf("\n");
	printf("\t\t\t-Tt:	Total simulation time in seconds.\n");
	printf("\t\t\t	[DEFAULT]  -Tt 1500\n");
	printf("\n");
	printf("\t\t\t-Tb:	Burn-in/threshold time. Time spent by the cell to reach equilibrium.\n");
	printf("\t\t\t	Only calculations after this time will be included in the analyses.\n");
	printf("\t\t\t	[DEFAULT]  -Tb 1000\n");
	printf("\n");
	printf("\t\t\t-Th:	Time at which harringtonine is added to the cell.\n");
	printf("\t\t\t	[DEFAULT]  -Th 1500\n");
	printf("\n");
	printf("\t\t\t-Tc:	Time at which cycloheximide is added to the cell.\n");
	printf("\t\t\t	[DEFAULT]  -Tc 1500\n");
	printf("\n");	
	printf("\t-R <value>	Total number of ribosomes in the cell.\n");
	printf("\t\t\t[DEFAULT]  -R 200000\n");
	printf("\n");
	printf("\t-t <value>	Total number of tRNAs in the cell.\n");
	printf("\t\t\t[DEFAULT]  -t 3300000\n");
	printf("\n");
	printf("\t-N <value>	Total number of genes. This needs to be specified by the user.\n");
	printf("\t\t\t[DEFAULT]  -N 1\n");
	printf("\n");
	printf("\t-F <FILE>	File containing processed fasta file into a numeric sequence.\n");
	printf("\t\t\tThis file is an output of the code utilities/convert.fasta.to.genom.pl\n");
	printf("\t\t\tIt contains the information regarding initiation probability, mRNA\n");
	printf("\t\t\tabundance and codon sequence of each gene.\n");
	printf("\t\t\t[DEFAULT]  -F example/input/S.cer.genom\n");
	printf("\n");
	printf("\t-C <FILE>	File containing the information about codon, tRNA, tRNA abundance and wobble.\n");
	printf("\t\t\tThis file is an output of the code utilities/create.Scer.cod.anticod.numeric.pl\n");
	printf("\t\t\t[DEFAULT]  -C example/input/S.cer.tRNA\n");
	printf("\n");
	printf("\t-J <FILE>	File containing the initial state of the system to begin simulations from.\n");
	printf("\t\t\tThis file is an output of this simulation code containing '*_ribo_pos_*'\n");
	printf("\t\t\t[DEFAULT]  -C example/input/output_final_ribo_pos.out\n");
	printf("\n");
	printf("\t-x[INTEGER]	Specify parameters for cycloheximide action\n");
	printf("\n");
	printf("\t\t\t-x1:	Probability of cycloheximide binding ribosomes\n");
	printf("\t\t\t	[DEFAULT]  -x1 0\n");
	printf("\n");
	printf("\t\t\t-x2:	Rate of cycloheximide dissociation from bound ribosomes\n");
	printf("\t\t\t	[DEFAULT]  -x2 0 (irreversible binding)\n");
	printf("\n");
	printf("\t-y <value>	Harringtonine rate for free ribosomes. This needs to be specified by the user.\n");
	printf("\t\t\t[DEFAULT]  -y 0\n");
	printf("\n");
	printf("\t-s <value>	Random number seed. *MUST SETUP*\n");
	printf("\t\t\t[DEFAULT]  -s 1\n");
	printf("\n");
	printf("\t-O <prefix>	Specifies the prefix for the output files.\n");
	printf("\t\t\t[DEFAULT]  -O output\n");
	printf("\n");
	printf("\n");
	printf("\t-p[INTEGER]	Specify which output files to print\n");
	printf("\n");
	printf("\t\t\t-p1:	Generates a file of average elongation times\n");
	printf("\t\t\t	of all codons.\n");
	printf("\n");
	printf("\t\t\t-p2:	Generates a file of total average elongation\n");
	printf("\t\t\t	time of each gene.\n");
	printf("\n");
	printf("\t\t\t-p3:	Generates a file of average time between initiation\n");
	printf("\t\t\t	events on mRNAs of each gene.\n");
	printf("\n");
	printf("\t\t\t-p4:	Generates a file of average number of free ribosomes,\n");
	printf("\t\t\t	and free tRNAs of each type.\n");
	printf("\n");
	printf("\t\t\t-p5:	Generates a file of the final state of all mRNAs in a cell.\n");
	printf("\t\t\t	It contains the poistions of all bound ribosomes on mRNAs.\n");
	printf("\n");
	printf("\t\t\t-p6:	This generates two files:\n");
	printf("\t\t\t	A file containing the amount of time wasted by stalled\n");
	printf("\t\t\t	ribosomes on mRNAs of each gene.\n");
	printf("\t\t\t	A file containing the time wasted by stalled ribosomes\n");
	printf("\t\t\t	on each codon position of Gene 0 (first gene in the\n");
	printf("\t\t\t	processed fasta file).\n");
	printf("\t\t\t	This option significantly increases the total running time\n");
	printf("\t\t\t	of the simulation. Use it with caution.\n");
	printf("\n");
	printf("\t\t\t-p7:	Generates a file of state of all mRNAs in a cell every second.\n");
	printf("\t\t\t	This is similar to -p5 printed every second.\n");
	printf("\t\t\t	This option significantly increases the total running time\n");
	printf("\t\t\t	of the simulation. Use it with caution.\n");
	printf("\n");
	printf("\t\t\t-p8:	Generates afile of number of bound ribosomes at each position of a gene.\n");
	printf("\n");
	printf("\t\t\t-p9:	Generates a file of average (RPF and mRNA based) of bound ribosomes\n");
	printf("\t\t\t	at each position of a gene.\n\n");
}

// Read in commandline arguments
void Read_Commandline_Args(int argc, char *argv[])
{	int i;
	
	for(i=1;i<argc;i++)
	{	if(argv[i][0] == '-')
		{	switch(argv[i][1])
			{	case 'h':
				case '-':
					fflush(stdout);
					Help_out();
					exit(1);
					break;
				case 'V':
					if((argv[i][2] != '\0') || (i==argc-1))
					{	printf("\nTotal space not specified or Incorrect usage\n");
						fflush(stdout);
						Help_out();
						exit(1);
					}
					else
					{	tot_space = atof(argv[++i]);
						
						avail_space_t = (double)floor(tot_space/pow(char_len_tRNA,3));
						avail_space_r = (double)floor(tot_space/pow(char_len_ribo,3));
						
						if(avail_space_r<1e3 || avail_space_t < 1e3)
						{	printf("\nAvailable cytoplasmic space is too small\n\n");
							fflush(stdout);
							Help_out();
							exit(1);
						}
						break;
					}
				case 's':
					if((argv[i][2] != '\0') || (i==argc-1))
					{	printf("\nSeed for RNG not specified or Incorrect usage\n");
						fflush(stdout);
						Help_out();
						exit(1);
					}
					else
					{	seed = atoi(argv[++i]);
						break;
					}
				case 't':
					if((argv[i][2] != '\0') || (i==argc-1))
					{	printf("\nTotal # of tRNAs not specified or Incorrect usage\n");
						fflush(stdout);
						Help_out();
						exit(1);
					}
					else
					{	tot_tRNA = atoi(argv[++i]);
						break;
					}
				case 'R':
					if((argv[i][2] != '\0') || (i==argc-1))
					{	printf("\nTotal # of ribosomes not specified or Incorrect usage\n");
						fflush(stdout);
						Help_out();
						exit(1);
					}
					else
					{	tot_ribo = atoi(argv[++i]);
						break;
					}
				case 'N':
					if((argv[i][2] != '\0') || (i==argc-1))
					{	printf("\nTotal # of genes not specified or Incorrect usage\n");
						fflush(stdout);
						Help_out();
						exit(1);
					}
					else
					{	n_genes = atoi(argv[++i]);
						if(n_genes>MAX_GENES)
						{	printf("\nTotal # of genes for simulation exceeds maximum genes = %d\n", MAX_GENES);
							fflush(stdout);
							Help_out();
							exit(1);
						}
						break;
					}
				case 'O':
					if((argv[i][2] != '\0') || (i==argc-1))
					{	printf("\nOutput prefix not specified or Incorrect usage\n");
						fflush(stdout);
						Help_out();
						exit(1);
					}
					else
					{	strcpy(out_prefix,argv[++i]);
						break;
					}
				case 'F':
					if((argv[i][2] != '\0') || (i==argc-1))
					{	printf("\nFasta file not specified or Incorrect usage\n");
						fflush(stdout);
						Help_out();
						exit(1);
					}
					else
					{	strcpy(fasta_file,argv[++i]);
						break;
					}
				case 'J':
					if((argv[i][2] != '\0') || (i==argc-1))
					{	printf("\nState file not specified or Incorrect usage\n");
						fflush(stdout);
						Help_out();
						exit(1);
					}
					else
					{	strcpy(state_file,argv[++i]);
						break;
					}
				case 'C':
					if((argv[i][2] != '\0') || (i==argc-1))
					{	printf("\nS.cer Code file not specified or Incorrect usage\n");
						fflush(stdout);
						Help_out();
						exit(1);
					}
					else
					{	strcpy(code_file,argv[++i]);
						break;
					}
				case 'x':
					switch(argv[i][2])
					{	case '1':
							chx_rb_prob = atof(argv[++i]);
							
							if(chx_rb_prob<0)
							{	printf("\nCycloheximide probability for bound ribosomes should be positive (>0)\n\n");
								fflush(stdout);
								Help_out();
								exit(1);
							}
							break;
						case '2':
							chx_ds_rate = atof(argv[++i]);
							
							if(chx_ds_rate<0)
							{	printf("\nCycloheximide dissociation rate from bound ribosomes should be positive (>0)\n\n");
								fflush(stdout);
								Help_out();
								exit(1);
							}
							break;
						default:
							printf("\nInvalid CHX options\n");
							fflush(stdout);
							Help_out();
							exit(1);
							break;
					}
					break;
				case 'y':
					if((argv[i][2] != '\0') || (i==argc-1))
					{	printf("\nHarringtonine rate for free ribosomes not specified or Incorrect usage\n");
						fflush(stdout);
						Help_out();
						exit(1);
					}
					else
					{	harr_rf_rate = atof(argv[++i]);
						
						if(harr_rf_rate<0)
						{	printf("\nHarringtonine rate for free ribosomes should be positive (>0)\n\n");
							fflush(stdout);
							Help_out();
							exit(1);
						}
						break;
					}
				case 'P':
				case 'p':
					switch(argv[i][2])
					{	case '1':
							printOpt[0]=1;		// Elongation times of all codons
							break;
						case '2':
							printOpt[1]=1;		// Average total elongation time of all genes
							break;
						case '3':
							printOpt[2]=1;		// Average time between initiation events of all genes
							break;
						case '4':
							printOpt[3]=1;		// Average number of free ribosomes and tRNAs of each type
							break;
						case '5':
							printOpt[4]=1;		// Final state of the system - positions on mRNAs bound by ribosomes
							break;
						case '6':
							printOpt[5]=1;		// Amount of time spent by ribosomes stalled on each gene
							break;
						case '7':
							printOpt[6]=1;		// Print internal state of the cell every second - positions on mRNAs bound by ribosomes
							break;
						case '8':
							printOpt[7]=1;		// Final state of the system - positions on genes bound by ribosomes
							break;
						case '9':
							printOpt[8]=1;		// Final state of the system - avg (RPF and mRNA normalized) ribosome densities
							break;
						default:
							printf("\nInvalid print options\n");
							fflush(stdout);
							Help_out();
							exit(1);
							break;
					}
				case 'T':
					switch(argv[i][2])
					{	case 't':
							tot_time = atof(argv[++i]);
							if(tot_time>MAX_TIME)
							{	printf("\nTotal time for simulation exceeds maximum allowed time = %d\n",MAX_TIME);
								fflush(stdout);
								Help_out();
								exit(1);
							}
							break;
						case 'b':
							thresh_time = atof(argv[++i]);
							if(thresh_time<0)
							{	printf("\nTime threshold %g should be > 0 and < Total time %g\n",thresh_time,tot_time);
								fflush(stdout);
								Help_out();
								exit(1);
							}
							break;
						case 'c':
							chx_time = atof(argv[++i]);
							if(chx_time<0)
							{	printf("\nTime for cycloheximide action %g should be > 0\n",chx_time);
								fflush(stdout);
								Help_out();
								exit(1);
							}
							break;
						case 'h':
							harr_time = atof(argv[++i]);
							if(harr_time<0)
							{	printf("\nTime for harringtonine action %g should be > 0\n",harr_time);
								fflush(stdout);
								Help_out();
								exit(1);
							}
							break;
					}
			}
		}
	}
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Main Function
int main(int argc, char *argv[])
{	int c1, c2, c3, c4;
	FILE *f1, *f2, *f3, *f4, *f5, *f6, *f7, *f8, *f9, *f10, *f11;

	// Read in arguments from the commandline
	Read_Commandline_Args(argc, argv);
	
	// Random number generation setup
	gsl_rng * r;
	gsl_rng_env_setup();
	
	r = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (r, (unsigned long) seed);

	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// User specified parameters for quick test
	gene *Gene = (gene *)malloc(sizeof(gene) * n_genes);
	if(Gene == NULL)
	{	printf("Too many genes\nOut of memory\n");fflush(stdout);
		exit(1);
	}
	
	// Initialize the various structures
	ribosome *Ribo = (ribosome *)malloc(sizeof(ribosome) * tot_ribo);
	if(Ribo == NULL)
	{	printf("Too many ribosomes\nOut of memory\n");fflush(stdout);
		exit(1);
	}
		
	trna *cTRNA = (trna *)malloc(sizeof(trna) * 61);
	if(cTRNA == NULL)
	{	printf("Too many tRNAs\nOut of memory\n");fflush(stdout);
		exit(1);
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Read in the numeric seq
	Read_FASTA_File(fasta_file, Gene);

	// Read in the trna code file
	Read_tRNA_File(code_file, cTRNA);
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	tot_mRNA = 0;
	for(c1=0;c1<n_genes;c1++)
	{	tot_mRNA += Gene[c1].exp;
	}
	
	transcript *mRNA = (transcript *)malloc(sizeof(transcript) * tot_mRNA);
	if(mRNA == NULL)
	{	printf("Too many mRNAs\nOut of memory\n");fflush(stdout);
		exit(1);
	}

	c3=0;
	for(c1=0;c1<n_genes;c1++)					// Fill the mRNA struct with
	{	for(c2=0;c2<Gene[c1].exp;c2++)			// length and ini_rate info.
		{	mRNA[c3].gene = c1;
			mRNA[c3].ini_n = 0;
			mRNA[c3].last_ini = 0.0;
			mRNA[c3].trans_n = 0;
			mRNA[c3].avg_time_to_ini = 0.0;
			mRNA[c3].avg_time_to_trans = 0.0;
			c3++;
		}
	}
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	// Variables to track the translation process
	int Tf[61];									// Number of free tRNAs of each type
	int Rf = tot_ribo;							// Total number of free ribosomes (at start all ribosomes are free)
	int Rfh = 0;								// Number of free ribosomes bound by harringtonine
	int Rb[61];									// Number of ribosomes bound to each codon
	int Mf[n_genes];							// Number of initiable mRNAs of each gene (first 10 codons unbound by ribosomes)
	int tot_Mf;									// Total number of initiable mRNAs
	int x;
	int n_Rb_e[61];								// Number of elongable bound ribosomes to each codon
	int r_id;
	int m_id;
	int c_id;
	int c2_id;
	int g_id;
	int tot_gcn=0;
	int obs_max_len = 0;						// Observed max gene length
	int obs_max_exp = 0;						// Observed max gene expression
	int next_avail_ribo = 0;
	int termtn_now=0;
	int n_e_times[61];							// Number of times a codon type is elongated
	double e_times[61];
	int num_waste_ribo[n_genes];				// Number of stalled ribosomes on mRNAs of each gene
	int num_waste_ribo_pos[Gene[0].len];		// Number of stalled ribosomes by codon position on mRNAs of the first gene
	double time_waste_ribo[n_genes];			// Total time spent by stalled ribosomes on mRNAs of each gene
	double time_waste_ribo_pos[Gene[0].len];	// Total time spent by stalled ribosomes at each codon position on mRNAs of the first gene
	double scld_Mf[n_genes];
	double tot_scld_Mf = 0.0;
	double r_ini;
	double r_elng[61];
	double avg_tRNA_abndc[61];					// Average number of free tRNAs of each type (averaged by time)
	double avg_Rf = {0.0};						// Average number of free ribosomes (averaged by time)
	
	int tot_stall_chx = 0;						// Total number of stalled ribosomes on mRNAs due to cycloheximide
	int num_stall_chx[n_genes];					// Number of stalled ribosomes on mRNAs of each gene due to cycloheximide
	int chx_r_id[tot_ribo];
	double cod_chx[61];							// Codon specific rate of cycloheximide action
	double tot_chx_ds_rate = 0.0;				// Total rate of CHX dissociaction
	double prob_chx_ds;

	double tot_harr_rate = 0.0;					// Total rate of harringtonine action
	double prob_harr;

	
	// Begin simulation of the translation process
	double t = 0.0;
	double prob_ini;
	double tmp_elng_prob;
	double prob_e[61];
	double prob_g;
	double tot_rate = 0.0;
	double inv_rate;
	double coin;
	int t_print = floor(thresh_time);
	char tmp_t[10];
	
	int **Rb_e;									// Bound ribosomes to each codon that can be elongated
	Rb_e = malloc(sizeof(int *) * 62);
	for(c1=0;c1<62;c1++)
	{	Rb_e[c1] = malloc(sizeof(int *) * tot_ribo);
	}
	
	int **gene_sp_ribo_500;
	gene_sp_ribo_500 = malloc(sizeof(int *) * n_genes);
	for(c1=0;c1<n_genes;c1++)
	{	gene_sp_ribo_500[c1] = malloc(sizeof(int *) * Gene[c1].len);
	}

	
	for(c1=0;c1<n_genes;c1++)
	{	if(Gene[c1].exp > obs_max_exp)
		{	obs_max_exp = Gene[c1].exp;
		}
		if(Gene[c1].len > obs_max_len)
		{	obs_max_len = Gene[c1].len;
		}
		Mf[c1] = Gene[c1].exp;
		scld_Mf[c1] = Mf[c1]*Gene[c1].ini_prob;
		tot_scld_Mf += scld_Mf[c1];							// Scale gene_exp with ini_prob
		
		num_waste_ribo[c1] = 0;								// Initialize stalled ribosomes and time spent
		num_stall_chx[c1] = 0;								// Initialize stalled ribosomes and time spent
		time_waste_ribo[c1] = 0.0;
	}

	for(c1=0;c1<Gene[0].len;c1++)							// Initialize stalled ribosomes and time spent at specific positions of gene 1
	{	num_waste_ribo_pos[c1] = 0;
		time_waste_ribo_pos[c1] = 0.0;
	}	
	
	int **free_mRNA;										// List to figure out which mRNAs can be initiated based on no bound ribosomes from pos=0->pos=10
	free_mRNA = malloc(sizeof(int *) * n_genes);
	for(c1=0;c1<n_genes;c1++)
	{	free_mRNA[c1] = malloc(sizeof(int *) * obs_max_exp);
	}

	// Initialize free mRNAs
	c3=0;
	for(c2=0;c2<n_genes;c2++)
	{	for(c1=0;c1<Gene[c2].exp;c1++)
		{	free_mRNA[c2][c1]=c3;							// At the begininning all mRNAs can be initiated (store their ids here)
			c3++;
		}
	}
	
	// Initialize free tRNAs and elongation time vectors
	for(c1=0;c1<61;c1++)
	{	if(cTRNA[c1].wobble==1.0)
		{	tot_gcn += cTRNA[c1].gcn;
		}
		e_times[c1] = 0.0;
		n_e_times[c1] = 0;
		
		Tf[c1]=0;
		avg_tRNA_abndc[c1]=0.0;
	}
	c2 = 0;
		
	for(c1=0;c1<61;c1++)
	{	if(cTRNA[c1].wobble==1.0)
		{	Tf[cTRNA[c1].tid] = floor((double)cTRNA[c1].gcn*tot_tRNA/tot_gcn);			// Number of tRNAs for tRNA-type cTRNA[c1].tid
		}
		n_Rb_e[c1] = 0;															// Initialize number of elongatable bound ribosomes to each codon
		cTRNA[c1].wobble = cTRNA[c1].wobble/(char_time_tRNA*avail_space_t);		// For faster computation
	}
	
	// Initialize R_grid
	// R_grid now contains the id of ribosome at each mRNA position.
	// If there is no ribosome then that position get the value tot_ribo instead of 0 as 0 is a ribosome id
	int **R_grid;							// The state of the system with respect to mRNAs and bound ribosomes
	R_grid = malloc(sizeof(int *) * tot_mRNA);
	for(c1=0;c1<tot_mRNA;c1++)
	{	R_grid[c1] = malloc(sizeof(int *) * obs_max_len);
		for(c2=0;c2<Gene[mRNA[c1].gene].len;c2++)
		{	R_grid[c1][c2] = tot_ribo;
		}
	}

	// If state file is provided, reinitialize the variables and arrays
	if(strcmp(state_file, "") != 0)												// Read in the state of the system
	{	Read_STATE_File(state_file, R_grid);

		m_id = 0;																// mRNA id
		r_id = 0;																// Ribosome id
		for(c1=0;c1<n_genes;c1++)
		{	for(c4=0;c4<Gene[c1].exp;c4++)
			{	for(c2=0;c2<10;c2++)											// For first 10 codons - check if mRNA is initiable
				{	if(R_grid[m_id][c2]==0)
					{	R_grid[m_id][c2] = tot_ribo;
					}
					else
					{	R_grid[m_id][c2] = r_id;
						
						Ribo[r_id].mRNA = m_id;
						Ribo[r_id].pos = c2;
						Ribo[r_id].inhbtr_bound = 0;
						
						c_id = Gene[mRNA[m_id].gene].seq[c2];					// Codon identity
						if(c2>0)
						{	Tf[cTRNA[Gene[mRNA[m_id].gene].seq[c2-1]].tid]--;
						}
						
						if(R_grid[m_id][c2+10]==0)								// Check if the current ribosome can be elongated
						{	Rb_e[c_id][n_Rb_e[c_id]] = r_id;					// If no ribosome at pos+10 then it can
							Ribo[r_id].elng_cod_list = c_id;
							Ribo[r_id].elng_pos_list = n_Rb_e[c_id];
							n_Rb_e[c_id]++;
						}
						else
						{	num_waste_ribo[mRNA[m_id].gene]++;
							if(mRNA[m_id].gene==0)
							{	num_waste_ribo_pos[c2]++;
							}
						}
						
						Mf[c1]--;
						scld_Mf[c1]-=Gene[c1].ini_prob;							// These mRNAs are not initiable
						tot_scld_Mf-=Gene[c1].ini_prob;
						if(c4!=Mf[c1])
						{	free_mRNA[c1][c4] = free_mRNA[c1][Mf[c1]];			// Accounting of free mRNAs for next round
						}
						
						r_id++;
						Rf--;
					}
				}
				for(c2=10;c2<Gene[mRNA[m_id].gene].len;c2++)
				{	if(R_grid[m_id][c2]==0)
					{	R_grid[m_id][c2] = tot_ribo;
					}
					else
					{	R_grid[m_id][c2] = r_id;
						
						Ribo[r_id].mRNA = m_id;
						Ribo[r_id].pos = c2;
						Ribo[r_id].inhbtr_bound = 0;
						
						c_id = Gene[mRNA[m_id].gene].seq[c2];									// Codon identity
						Tf[cTRNA[Gene[mRNA[m_id].gene].seq[c2-1]].tid]--;
						
						if(R_grid[m_id][c2+10]==0 || (c2+11)>=Gene[mRNA[m_id].gene].len)		// Check if the current ribosome can be elongated
						{	Rb_e[c_id][n_Rb_e[c_id]] = r_id;									// If no ribosome at pos+10 then it can
							Ribo[r_id].elng_cod_list = c_id;
							Ribo[r_id].elng_pos_list = n_Rb_e[c_id];
							n_Rb_e[c_id]++;
						}
						else
						{	num_waste_ribo[mRNA[m_id].gene]++;
							if(mRNA[m_id].gene==0)
							{	num_waste_ribo_pos[c2]++;
							}
						}
						
						r_id++;
						Rf--;
					}
				}
				m_id++;
			}
		}
		
		next_avail_ribo = r_id;
	}
	

	/////////////////////////////////////////////////
	// Begin the actual continuous time simulation
	/////////////////////////////////////////////////
	
	while(t<tot_time)																// Till current time is less than max time
	{	r_ini = tot_scld_Mf*Rf/(char_time_ribo*avail_space_r);						// Initiation rate
		tot_rate = r_ini;
		
		for(c1=0;c1<61;c1++)
		{	r_elng[c1] = Tf[cTRNA[c1].tid]*cTRNA[c1].wobble;						// Elongation rate of codon c1
			r_elng[c1] *= (double)n_Rb_e[c1];
			tot_rate += r_elng[c1];
		}
		
		if(t>harr_time)
		{	tot_harr_rate = (double)Rf*harr_rf_rate;								// Rate of harringtonine action on free ribosomes
			tot_rate += tot_harr_rate;
		}

		if(t>chx_time)
		{	tot_chx_ds_rate = (double)tot_stall_chx*chx_ds_rate;					// CHX dissociation rate action from bound ribosomes
			tot_rate += tot_chx_ds_rate;
		}

		if(tot_rate>0)
		{	inv_rate = 1/tot_rate;
		}
		else
		{	printf("\nNo further events to process.\nSimulation stopped at time %g\n\n",t);fflush(stdout);
			t=tot_time;
			continue;
		}
		
		// Increment time
		if(t>thresh_time)
		{	if(printOpt[3]==1)														// Tracking free ribosomes and tRNAs
			{	for(c1=0;c1<61;c1++)
				{	avg_tRNA_abndc[c1] += (double)Tf[c1]*inv_rate;
				}
				avg_Rf += (double)Rf*inv_rate;
			}
			
			if(printOpt[5]==1)														// Tracking time wasted by ribosome stalling
			{	for(c1=0;c1<n_genes;c1++)
				{	time_waste_ribo[c1] += (double)num_waste_ribo[c1]*inv_rate;
				}
				for(c1=0;c1<Gene[0].len;c1++)
				{	time_waste_ribo_pos[c1] += (double)num_waste_ribo_pos[c1]*inv_rate;
				}
			}
		}
		t+=inv_rate;
		
		// Calculate prob of events
		prob_ini = r_ini*inv_rate;													// Prob of initiation in current t
		prob_harr = prob_ini + tot_harr_rate*inv_rate;								// Prob of harr action in current t
		prob_chx_ds = prob_harr + tot_chx_ds_rate*inv_rate;							// Prob of CHX dissociaction in current t
			
		coin = gsl_rng_uniform(r);													// Pick a random uniform to pick an event

		// Translation initiation
		if(coin<prob_ini)
		{	
			// Ribosomes are picked sequentially
			r_id = next_avail_ribo;
			Ribo[r_id].pos = 0;
			Ribo[r_id].t_trans_ini = t;
			Ribo[r_id].t_elong_ini = t;
			next_avail_ribo++;
			
			// Pick a random mRNA for initiation
			coin = gsl_rng_uniform(r);
			prob_g = 0.0;
			
			for(c1=0;c1<n_genes;c1++)
			{	prob_g += scld_Mf[c1]/tot_scld_Mf;									// Pick a gene randomly first as they may differ in ini_prob
				if(coin<prob_g)
				{	c2 = gsl_rng_uniform_int(r, (unsigned long)Mf[c1]);				// Once a gene is selected pick a random mRNA
					m_id = free_mRNA[c1][c2];
					Mf[c1]--;
					scld_Mf[c1]-=Gene[c1].ini_prob;
					tot_scld_Mf-=Gene[c1].ini_prob;
					if(c2!=Mf[c1])
					{	free_mRNA[c1][c2] = free_mRNA[c1][Mf[c1]];					// Accounting of free mRNAs for next round
					}
					break;
				}
			}			
			
			if(t>thresh_time)
			{	mRNA[m_id].ini_n++;													// Store the # of initn events on this mRNA
				mRNA[m_id].avg_time_to_ini += t-mRNA[m_id].last_ini;				// Time between initn
			}	
			mRNA[m_id].last_ini = t;

			c_id = Gene[mRNA[m_id].gene].seq[0];									// Codon at pos 1

			Ribo[r_id].mRNA = m_id;							
			
			coin = gsl_rng_uniform(r);
			prob_g = (double)Rfh/(double)Rf;
			if(coin>prob_g)
			{	if(R_grid[m_id][10]==tot_ribo)										// Check if the current ribosome can be elongated
				{	Rb_e[c_id][n_Rb_e[c_id]] = r_id;								// If no ribosome at pos+10 then it can
					Ribo[r_id].elng_cod_list = c_id;
					Ribo[r_id].elng_pos_list = n_Rb_e[c_id];
					n_Rb_e[c_id]++;
				}
				else
				{	num_waste_ribo[mRNA[m_id].gene]++;
					if(mRNA[m_id].gene==0)
					{	num_waste_ribo_pos[0]++;
					}
				}
				Ribo[r_id].inhbtr_bound = 0;
			}
			else
			{	num_waste_ribo[mRNA[m_id].gene]++;
				Rfh--;
				Ribo[r_id].inhbtr_bound = 1;
			}
			R_grid[m_id][0] = r_id;													// Update the ribosome grid uypon initiation
			
			Rf--;																	// Update number of free ribosomes
		}
		// Harringtonine action
		else if(coin<prob_harr)
		{	Rfh++;
		}
		// CHX dissociaction
		else if(coin<prob_chx_ds)
		{	x = gsl_rng_uniform_int(r, (unsigned long)tot_stall_chx);				// Randomly pick a CHX bound ribosome
			r_id = chx_r_id[x];
			m_id = Ribo[r_id].mRNA;
			c_id = Gene[mRNA[m_id].gene].seq[Ribo[r_id].pos];
			
			if((Ribo[r_id].pos+11)>=Gene[mRNA[m_id].gene].len || R_grid[m_id][Ribo[r_id].pos+11]==tot_ribo)
			{	Rb_e[c_id][n_Rb_e[c_id]] = r_id;
				Ribo[r_id].elng_cod_list = c_id;
				Ribo[r_id].elng_pos_list = n_Rb_e[c_id];
				n_Rb_e[c_id]++;
			}
			if(Ribo[r_id].pos!=(Gene[mRNA[m_id].gene].len-1) && Ribo[r_id].pos>1)	// Check if the current elongation has led to termination
			{	Tf[cTRNA[Gene[mRNA[m_id].gene].seq[Ribo[r_id].pos-2]].tid]++;
			}
			Ribo[r_id].inhbtr_bound = 0;
			
			tot_stall_chx--;
			chx_r_id[x] = chx_r_id[tot_stall_chx];
			num_stall_chx[mRNA[m_id].gene]--;
		}
		// Translation Elongation
		else
		{	tmp_elng_prob = prob_chx_ds;
			for(c1=0;c1<61;c1++)	
			{	prob_e[c1] = tmp_elng_prob + r_elng[c1]*inv_rate;					// Prob of elongn of codon 1
				if(coin<prob_e[c1])
				{	c_id = c1;
					break;
				}
				tmp_elng_prob = prob_e[c1];
			}
			
			x = gsl_rng_uniform_int(r, (unsigned long)n_Rb_e[c_id]);				// Randomly pick an elongatable ribosome bound to codon c_id
			r_id = Rb_e[c_id][x];
			m_id = Ribo[r_id].mRNA;
			
			// Elongation arrest due to Cycloheximide action
			coin = gsl_rng_uniform(r);
			if(coin<chx_rb_prob && t>chx_time)
			{	n_Rb_e[c_id]--;														// Remove it from the list of elongatable codons
				if(x!=n_Rb_e[c_id])
				{	Ribo[Rb_e[c_id][n_Rb_e[c_id]]].elng_pos_list = x;				// Update the ids and number of elongatable ribosomes
					Rb_e[c_id][x] = Rb_e[c_id][n_Rb_e[c_id]];
				}
				if(Ribo[r_id].pos!=(Gene[mRNA[m_id].gene].len-1))					// Check if the current elongation has led to termination
				{	Tf[cTRNA[c_id].tid]--;
					R_grid[m_id][Ribo[r_id].pos] = tot_ribo;
					Ribo[r_id].pos++;
					R_grid[m_id][Ribo[r_id].pos] = r_id;
				}
				else
				{	termtn_now=1;
				}
				num_stall_chx[mRNA[m_id].gene]++;

				Ribo[r_id].inhbtr_bound = 1;

				chx_r_id[tot_stall_chx] = r_id;
				tot_stall_chx++;
			}
			// Normal elongation cycle
			else
			{	if(Ribo[r_id].pos>0)													// Release the tRNA bound at the earlier position
				{	Tf[cTRNA[Gene[mRNA[m_id].gene].seq[Ribo[r_id].pos-1]].tid]++;
				}
						
				if(Ribo[r_id].pos==(Gene[mRNA[m_id].gene].len-1))						// Check if the current elongation has led to termination
				{	R_grid[m_id][Ribo[r_id].pos] = tot_ribo;							// Update the position of ribosomes on the mRNA
				
					Ribo[r_id].pos++;
					Rf++;																// Free a ribosome upon termination

					n_Rb_e[c_id]--;
					if(x!=n_Rb_e[c_id])
					{	Ribo[Rb_e[c_id][n_Rb_e[c_id]]].elng_pos_list = x;				// Update the ids and number of elongatable ribosomes
						Rb_e[c_id][x] = Rb_e[c_id][n_Rb_e[c_id]];
					}
			
					if(t>thresh_time)
					{	if(printOpt[0]==1)
						{	e_times[c_id] += t-Ribo[r_id].t_elong_ini;					// For estimation of avg elongation times of codons
							n_e_times[c_id]++;
						}
						mRNA[m_id].trans_n++;											// Update the number of trans evnts on curr mRNA
						mRNA[m_id].avg_time_to_trans += t-Ribo[r_id].t_trans_ini;		// Update the time to translation
					}

					// Update any previously unelongatable ribosomes
					if(R_grid[m_id][Ribo[r_id].pos-11]<tot_ribo)						// When the ribosome moves, a previously unelongatable ribosome
					{	c2_id = Gene[mRNA[m_id].gene].seq[Ribo[r_id].pos-11];			// can now be elongatable on the same mRNA if its 11 codon behind
						Rb_e[c2_id][n_Rb_e[c2_id]] = R_grid[m_id][Ribo[r_id].pos-11];
					
						Ribo[R_grid[m_id][Ribo[r_id].pos-11]].elng_cod_list = c2_id;
						Ribo[R_grid[m_id][Ribo[r_id].pos-11]].elng_pos_list = n_Rb_e[c2_id];
						n_Rb_e[c2_id]++;

						num_waste_ribo[mRNA[m_id].gene]--;
						if(mRNA[m_id].gene==0)
						{	num_waste_ribo_pos[Ribo[r_id].pos-11]--;
						}
					}
								
					// Update the ribosomes avail for initiation
					next_avail_ribo--;
				
					if(r_id!=next_avail_ribo)											// The last initiated ribosome's data is swapped with
					{	Ribo[r_id].mRNA = Ribo[next_avail_ribo].mRNA;					// the ribosome that just finished translation (terminated)
						Ribo[r_id].pos = Ribo[next_avail_ribo].pos;						// This is done to primarily keep track of ONLY ribosomes
						Ribo[r_id].t_trans_ini = Ribo[next_avail_ribo].t_trans_ini;		// that are currently bound for faster computation
						Ribo[r_id].t_elong_ini = Ribo[next_avail_ribo].t_elong_ini;
						
						Ribo[r_id].elng_cod_list = Gene[mRNA[Ribo[next_avail_ribo].mRNA].gene].seq[Ribo[next_avail_ribo].pos];
						Ribo[r_id].elng_pos_list = Ribo[next_avail_ribo].elng_pos_list;
						R_grid[Ribo[r_id].mRNA][Ribo[r_id].pos] = r_id;
						
						if(R_grid[Ribo[r_id].mRNA][Ribo[r_id].pos+10]==tot_ribo || (Ribo[r_id].pos+10)>=Gene[mRNA[Ribo[r_id].mRNA].gene].len)
						{	Rb_e[Ribo[r_id].elng_cod_list][Ribo[r_id].elng_pos_list] = r_id;
						}
					}
					termtn_now = 1;
				
				}
				else if((Ribo[r_id].pos+11)>=Gene[mRNA[m_id].gene].len || R_grid[m_id][Ribo[r_id].pos+11]==tot_ribo)	// Check if the ribosome is still elongatable
				{	R_grid[m_id][Ribo[r_id].pos] = tot_ribo;															// Update the position of ribosomes on the mRNA
					
					Ribo[r_id].pos++;
					R_grid[m_id][Ribo[r_id].pos] = r_id;
					c2_id = Gene[mRNA[m_id].gene].seq[Ribo[r_id].pos];
					
					if(c2_id!=c_id)															// If the codon has changed shift the elongatable ribosome
					{	Rb_e[c2_id][n_Rb_e[c2_id]] = r_id;									// to the other codon
						Ribo[r_id].elng_cod_list = c2_id;
						Ribo[r_id].elng_pos_list = n_Rb_e[c2_id];
						n_Rb_e[c2_id]++;
						
						n_Rb_e[c_id]--;
						if(x!=n_Rb_e[c_id])
						{	Ribo[Rb_e[c_id][n_Rb_e[c_id]]].elng_pos_list = x;				// Update the ids and number of elongatable ribosomes
							Rb_e[c_id][x] = Rb_e[c_id][n_Rb_e[c_id]];
						}
					}
					Tf[cTRNA[c_id].tid]--;
				}
				else																		// When ribosome is not elongatable anymore
				{	R_grid[m_id][Ribo[r_id].pos] = tot_ribo;								// Update the position of ribosomes on the mRNA
				
					Ribo[r_id].pos++;
					R_grid[m_id][Ribo[r_id].pos] = r_id;
					Ribo[r_id].elng_cod_list = Gene[mRNA[m_id].gene].seq[Ribo[r_id].pos];

					n_Rb_e[c_id]--;
					if(x!=n_Rb_e[c_id])
					{	Ribo[Rb_e[c_id][n_Rb_e[c_id]]].elng_pos_list = x;					// Update the ids and number of elongatable ribosomes
						Rb_e[c_id][x] = Rb_e[c_id][n_Rb_e[c_id]];
					}
				
					Tf[cTRNA[c_id].tid]--;
					num_waste_ribo[mRNA[m_id].gene]++;
					if(mRNA[m_id].gene==0)
					{	num_waste_ribo_pos[Ribo[r_id].pos]++;
					}
				}
			}
			if(termtn_now==0)															// When elongation does not lead to termination
			{	if(Ribo[r_id].pos>10 && R_grid[m_id][Ribo[r_id].pos-11]<tot_ribo)		// When the ribosome moves, a previously unelongatable ribosome
				{	if(Ribo[R_grid[m_id][Ribo[r_id].pos-11]].inhbtr_bound==0)			// can now be elongated on the same mRNA if its 11 codon behind
					{	c2_id = Gene[mRNA[m_id].gene].seq[Ribo[r_id].pos-11];
						Rb_e[c2_id][n_Rb_e[c2_id]] = R_grid[m_id][Ribo[r_id].pos-11];
					
						Ribo[R_grid[m_id][Ribo[r_id].pos-11]].elng_cod_list = c2_id;
						Ribo[R_grid[m_id][Ribo[r_id].pos-11]].elng_pos_list = n_Rb_e[c2_id];
						n_Rb_e[c2_id]++;
					
						num_waste_ribo[mRNA[m_id].gene]--;
						if(mRNA[m_id].gene==0)
						{	num_waste_ribo_pos[Ribo[r_id].pos-11]--;
						}
					}
				}
				if(Ribo[r_id].pos==10)													// If the current position is at codon 11 make the current mRNA initiable
				{	g_id = mRNA[m_id].gene;
					free_mRNA[g_id][Mf[g_id]] = m_id;
					Mf[g_id]++;
					scld_Mf[g_id] += Gene[g_id].ini_prob;
					tot_scld_Mf += Gene[g_id].ini_prob;
				}

				Ribo[r_id].t_elong_ini = t;												// Upon elongation, update the elong ini time for the next evnt
			}
			else
			{	termtn_now = 0;
			}
			
			if(t>((double)t_print) && printOpt[6]==1)
			{	for(c1=0;c1<tot_mRNA;c1++)
				{	for(c2=0;c2<Gene[mRNA[c1].gene].len;c2++)
					{	if(R_grid[c1][c2]!=tot_ribo)
						{	gene_sp_ribo_500[mRNA[c1].gene][c2]++;
						}
					}
				}
				t_print++;
			}
		}
	}
	
	
	
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Process output for printing
		
	int n_trans[n_genes];
	int **gene_sp_ribo;
	gene_sp_ribo = malloc(sizeof(int *) * n_genes);
	for(c1=0;c1<n_genes;c1++)
	{	gene_sp_ribo[c1] = malloc(sizeof(int *) * Gene[c1].len);
	}
	int **gene_sp_ribo_ns;
	gene_sp_ribo_ns = malloc(sizeof(int *) * n_genes);
	for(c1=0;c1<n_genes;c1++)
	{	gene_sp_ribo_ns[c1] = malloc(sizeof(int *) * Gene[c1].len);
	}
	
	int n_tot=0;
	int n_tot_ns=0;
	int n_tot_all=0;
	int gcount[obs_max_len];
	double gmean = 0.0;
	double gmean_ns = 0.0;
	double avg_ribo_pos[4][obs_max_len];
	double g_etimes[n_genes];
	double g_ini[n_genes];
	double dCST=0.0;
	
	if(printOpt[1]==1 || printOpt[2]==1 || printOpt[3]==1)
	{	c3 = 0;
		for(c1=0;c1<n_genes;c1++)
		{	n_trans[c1] = 0;
			g_ini[c1]=0.0;
			g_etimes[c1]=0.0;
			
			for(c2=0;c2<Gene[c1].exp;c2++)
			{	n_trans[c1] += mRNA[c3].trans_n;
				g_ini[c1] += mRNA[c3].avg_time_to_ini;
				g_etimes[c1] += mRNA[c3].avg_time_to_trans;
				n_tot += mRNA[c3].trans_n;
				c3++;
			}
		}
	}
	
	
	// Print the output
	
	// Elongation times of all codons
	if(printOpt[0]==1)
	{	strcpy(out_file,out_prefix);
		f2 = fopen(strcat(out_file,"_etimes.out"),"w");
		
		fprintf(f2,"Codon\tNum_of_events\tAvg_elong_time(sec)\n");
		for(c1=0;c1<61;c1++)
		{	dCST = e_times[c1]/(double)n_e_times[c1];
			fprintf(f2,"%d\t%d\t%g\n",c1,n_e_times[c1],dCST);
		}
		fclose(f2);
	}
	
	// Average total elongation times of all genes
	if(printOpt[1]==1)
	{	strcpy(out_file,out_prefix);
		f3 = fopen(strcat(out_file,"_gene_totetimes.out"),"w");
		
		fprintf(f3,"Gene\tNum_of_events\tAvg_total_elong_time(sec)\tNum_elng_stall\tNum_chx_stall\n");
		for(c1=0;c1<n_genes;c1++)
		{	dCST = g_etimes[c1]/(double)n_trans[c1];
			fprintf(f3,"%d\t%d\t%g\t%d\t%d\n",c1,n_trans[c1],dCST,num_waste_ribo[c1],num_stall_chx[c1]);
		}
		fclose(f3);
	}
	
	// Average time between initiation of all genes
	if(printOpt[2]==1)
	{	strcpy(out_file,out_prefix);
		f4 = fopen(strcat(out_file,"_gene_initimes.out"),"w");
		
		fprintf(f4,"Gene\tNum_of_events\tAvg_initiation_time(sec)\n");
		for(c1=0;c1<n_genes;c1++)
		{	dCST = g_ini[c1]/(double)n_trans[c1];
			fprintf(f4,"%d\t%d\t%g\n",c1,n_trans[c1],dCST);
		}
		fclose(f4);
	}
	
	// Average number of free ribosomes and tRNAs at equilibrium
	if(printOpt[3]==1)
	{	strcpy(out_file,out_prefix);
		f5 = fopen(strcat(out_file,"_avg_ribo_tRNA.out"),"w");
		
		avg_Rf = avg_Rf/(tot_time-thresh_time);
		fprintf(f5,"Free_ribo\t%g\n",avg_Rf);
		for(c1=0;c1<61;c1++)
		{	if(avg_tRNA_abndc[c1]>0)
		{	avg_tRNA_abndc[c1] = avg_tRNA_abndc[c1]/(tot_time-thresh_time);
			fprintf(f5,"Free_tRNA%d\t%g\n",c1,avg_tRNA_abndc[c1]);
		}
		}
		fclose(f5);
	}
	
	// The final state of the system - positions of bound ribosomes on mRNAs
	if(printOpt[4]==1)
	{	strcpy(out_file,out_prefix);
		f6 = fopen(strcat(out_file,"_final_ribo_pos.out"),"w");
		
		// Print final state for individual mRNAs
		for(c1=0;c1<tot_mRNA;c1++)
		{	if(R_grid[c1][0]==tot_ribo)
		{	fprintf(f6,"0");
		}
		else
		{	fprintf(f6,"1");
		}
			
			for(c2=1;c2<Gene[mRNA[c1].gene].len;c2++)
			{	if(R_grid[c1][c2]==tot_ribo)
			{	fprintf(f6," 0");
			}
			else
			{	fprintf(f6," 1");
			}
			}
			fprintf(f6,"\n");
		}
		fclose(f6);
	}
	
	// Time spent by stalled ribosomes on each gene
	if(printOpt[5]==1)
	{	strcpy(out_file,out_prefix);
		f7 = fopen(strcat(out_file,"_gene0_pos_stall_ribo.out"),"w");
		
		fprintf(f7,"Pos\tAvg_ribo_stall\n");
		for(c1=0;c1<Gene[0].len;c1++)
		{	time_waste_ribo_pos[c1] = time_waste_ribo_pos[c1]/(tot_time-thresh_time);
			fprintf(f7,"%d\t%g\n",c1,time_waste_ribo_pos[c1]);
		}
		fclose(f7);
		
		strcpy(out_file,out_prefix);
		f8 = fopen(strcat(out_file,"_allgene_stall_ribo.out"),"w");
		fprintf(f8,"Gene\tAvg_ribo_stall\n");
		for(c1=0;c1<n_genes;c1++)
		{	time_waste_ribo[c1] = time_waste_ribo[c1]/(tot_time-thresh_time);
			fprintf(f8,"%d\t%g\n",c1,time_waste_ribo[c1]);
		}
		fclose(f8);
	}
	
	
	if(printOpt[6]==1)
	{	strcpy(out_file,out_prefix);
		f9 = fopen(strcat(out_file,"_gene_pos_sp_ribo_500.out"),"w");
		
		for(c1=0;c1<n_genes;c1++)
		{	fprintf(f9,"%d",gene_sp_ribo_500[c1][0]);
			for(c2=1;c2<Gene[c1].len;c2++)
			{	fprintf(f9," %d",gene_sp_ribo_500[c1][c2]);
			}
			fprintf(f9,"\n");
		}
		fclose(f9);
	}
	
	
	// The final state of the system - positions of bound ribosomes on mRNAs and genes
	if(printOpt[7]==1 || printOpt[8]==1)
	{	// Initialize the gene position specific ribosome density array
		for(c1=0;c1<n_genes;c1++)
		{	for(c2=0;c2<Gene[c1].len;c2++)
			{	gene_sp_ribo[c1][c2] = 0;
				gene_sp_ribo_ns[c1][c2] = 0;
			}
		}
		
		// Initialize the avg position specific ribosome density array
		for(c1=0;c1<obs_max_len;c1++)
		{	avg_ribo_pos[0][c1] = 0.0;
			avg_ribo_pos[1][c1] = 0.0;
			avg_ribo_pos[2][c1] = 0.0;
			avg_ribo_pos[3][c1] = 0.0;
			gcount[c1] = 0;
		}
		
		// Final state for individual mRNAs
		for(c1=0;c1<tot_mRNA;c1++)
		{	if(R_grid[c1][0]!=tot_ribo)
			{	gene_sp_ribo[mRNA[c1].gene][0]++;
				if(R_grid[c1][10]==tot_ribo)
				{	gene_sp_ribo_ns[mRNA[c1].gene][0]++;
				}
			}
			
			for(c2=1;c2<Gene[mRNA[c1].gene].len;c2++)
			{	if(R_grid[c1][c2]!=tot_ribo)
				{	gene_sp_ribo[mRNA[c1].gene][c2]++;
					if(R_grid[c1][c2+10]==tot_ribo || (c2+10)>Gene[mRNA[c1].gene].len)
					{	if(R_grid[c1][c2-10]==tot_ribo || c2<10)
						{	gene_sp_ribo_ns[mRNA[c1].gene][c2]++;
						}
					}
				}
			}
		}
		
		// Calculate final density on each position of a gene
		for(c1=0;c1<n_genes;c1++)
		{	n_tot = 0;
			n_tot_ns = 0;
			
			gcount[0]++;
			n_tot += gene_sp_ribo[c1][0];
			n_tot_ns += gene_sp_ribo_ns[c1][0];
			
			for(c2=1;c2<Gene[c1].len;c2++)
			{	gcount[c2]++;
				n_tot += gene_sp_ribo[c1][c2];
				n_tot_ns += gene_sp_ribo_ns[c1][c2];
			}
			
			gmean = (double)n_tot/(double)Gene[c1].len;
			gmean_ns = (double)n_tot_ns/(double)Gene[c1].len;
			for(c2=0;c2<Gene[c1].len;c2++)
			{	if(gmean>0)
				{	avg_ribo_pos[0][c2] += (double)gene_sp_ribo[c1][c2]/gmean;
				}
				if(gmean_ns>0)
				{	avg_ribo_pos[2][c2] += (double)gene_sp_ribo_ns[c1][c2]/gmean_ns;
				}
				avg_ribo_pos[1][c2] += (double)gene_sp_ribo[c1][c2]/(double)Gene[c1].exp;
				avg_ribo_pos[3][c2] += (double)gene_sp_ribo_ns[c1][c2]/(double)Gene[c1].exp;
			}
		}
		
		if(printOpt[7]==1)
		{	strcpy(out_file,out_prefix);
			f9 = fopen(strcat(out_file,"_gene_pos_sp_ribo.out"),"w");
			
			for(c1=0;c1<n_genes;c1++)
			{	fprintf(f9,"%d",gene_sp_ribo[c1][0]);
				for(c2=1;c2<Gene[c1].len;c2++)
				{	fprintf(f9," %d",gene_sp_ribo[c1][c2]);
				}
				fprintf(f9,"\n");
			}
			fclose(f9);
			
			strcpy(out_file,out_prefix);
			f11 = fopen(strcat(out_file,"_gene_pos_sp_ribo_nostall.out"),"w");
			
			for(c1=0;c1<n_genes;c1++)
			{	fprintf(f11,"%d",gene_sp_ribo_ns[c1][0]);
				for(c2=1;c2<Gene[c1].len;c2++)
				{	fprintf(f11," %d",gene_sp_ribo_ns[c1][c2]);
				}
				fprintf(f11,"\n");
			}
			fclose(f11);
		}
		
		if(printOpt[8]==1)
		{	strcpy(out_file,out_prefix);
			f10 = fopen(strcat(out_file,"_avg_pos_sp_ribo.out"),"w");
			fprintf(f10,"Pos\tRPF_avg\tmRNA_avg\tRPF_avg_nostall\tmRNA_avg_nostall\n");
			for(c2=0;c2<obs_max_len;c2++)
			{	c1 = c2+1;
				avg_ribo_pos[0][c2] /= (double)gcount[c2];
				avg_ribo_pos[1][c2] /= (double)gcount[c2];
				avg_ribo_pos[2][c2] /= (double)gcount[c2];
				avg_ribo_pos[3][c2] /= (double)gcount[c2];
				fprintf(f10,"%d\t%g\t%g\t%g\t%g\n",c1,avg_ribo_pos[0][c2],avg_ribo_pos[1][c2],avg_ribo_pos[2][c2],avg_ribo_pos[3][c2]);
			}
			fclose(f10);
		}
	}
	
	
	// Free the malloc structures and arrays
	free(Gene);
	free(Ribo);
	free(mRNA);
	free(cTRNA);
	free(Rb_e);
	free(gene_sp_ribo);
	free(gene_sp_ribo_ns);
	free(gene_sp_ribo_500);
	free(free_mRNA);
	free(R_grid);
	
	gsl_rng_free(r);
}
