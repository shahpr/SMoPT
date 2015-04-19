# Convert a fasta sequence to a numeric sequence based on codon numbers.
# This is for algorithmic convenience in both C and R, and leads to faster calculations.

# Usage: 	perl convert.fasta.to.genom.pl <fasta_file> <genetic_code_file> <mRNA_abndc_initiation_file> <output_file>
# Example:	perl convert.fasta.to.genom.pl ../example/input/S.cer.fasta genetic.code.tsv ../example/input/S.cer.mRNA.abndc.ini.tsv ../example/input/S.cer

# Read in fasta file
open fi,"$ARGV[0]";
chomp(@fas_file=<fi>);

# Read in genetic code file
open fi2,"$ARGV[1]";
chomp(@t=<fi2>);

# Read in file with list of genes, their initiation probabilities and mRNA abundances
open fi3,"$ARGV[2]";
chomp(@mR=<fi3>);

# Output for processed fasta file
open fo,">$ARGV[3].genom";

# Process the genetic code file
%trna=();
$c=0;
for($i=1;$i<@t;$i++)
{	@a=split(/\t/,$t[$i]);
	$trna{$a[1]}=$c;
	$c++;
}

# Process mRNA abundance and initiation prob file
@gpar=();
for($i=1;$i<@mR;$i++)
{	@a=split(/\s+/,$mR[$i]);
	$gpar[$i-1][0]=$a[1];
	$gpar[$i-1][1]=int($a[2]);
}

# Process and print fasta file
$c=0;
for($i=0;$i<@fas_file;$i++)
{	if($fas_file[$i]=~">")
	{	$i++;
		$tmp="";
        	do  
                {	$tmp.=$fas_file[$i];
			$i++;
		}while($fas_file[$i]!~">" & $i<@fas_file);
		$i--;
        }

	$l=length($tmp);

	print fo "$gpar[$c][0] $gpar[$c][1] ";
	for($j=0;$j<($l-6);$j+=3)
	{	$cod=substr($tmp,$j,3);
		print fo "$trna{$cod} ";
	}
	$cod=substr($tmp,$j,3);
	print fo "$trna{$cod}\n";

	$c++;
}
