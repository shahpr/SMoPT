# Prepare a tRNA file for simulation that contains the following information:
# tRNA type regonizing a given codon gene copy number of that tRNA and 
# the wobble parameter between the focal codon and the tRNA type.

# Read in file containing the mapping between codons and anticodons
open fi, "cod.anticod.file.tsv";
chomp(@s=<fi>);

# List of gene copy numbers corresponding to each tRNA type
open fi2,"S.cer.tRNA_copy_num";
chomp(@t=<fi2>);

# Mapping between amino acid and codons. (Serine is split into two sets)
open fi3,"genetic.code.tsv";
chomp(@f=<fi3>);

# Output tRNA file for simulation
open fo,">../example/input/S.cer.tRNA";

%gcn=();
$c=0;
for($i=1;$i<@t;$i++)
{	@a=split(/\t/,$t[$i]);
	if($a[1]>0)
	{	$gcn{$a[0]}=[$a[1],$c];
		$c++;
	}
}

%rln=();
foreach(@s)
{	@a=split(/\t/,$_);
	$rln{$a[0]}=[$a[1],$a[3]];
}

for($i=1;$i<@f;$i++)
{	@a=split(/\t/,$f[$i]);
	$j=$i-1;
	if($gcn{$rln{$a[1]}[0]}[0])													# If perfect match between codon and tRNA (wobble=1)
	{	print fo "$j\t$gcn{$rln{$a[1]}[0]}[1]\t$gcn{$rln{$a[1]}[0]}[0]\t1\n";
	}
	else
	{	$w=substr($a[1],2,1).substr($rln{$a[1]}[1],0,1);
		if($w eq "AA" || $w eq "AG" || $w eq "GA" || $w eq "GG" || $w eq "TT" || $w eq "TC" || $w eq "CT" || $w eq "CC")	# Wobble RR/YY = 0.61
		{	print fo "$j\t$gcn{$rln{$a[1]}[1]}[1]\t$gcn{$rln{$a[1]}[1]}[0]\t0.61\n";
		}
		else
		{	print fo "$j\t$gcn{$rln{$a[1]}[1]}[1]\t$gcn{$rln{$a[1]}[1]}[0]\t0.64\n";					# Wobble RY/YR = 0.64
		}
	}
}
