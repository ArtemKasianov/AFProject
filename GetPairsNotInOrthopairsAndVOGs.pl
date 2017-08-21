use strict;

my $identFile = $ARGV[0];
my $orthopairsFile = $ARGV[1];
my $vogsFile = $ARGV[2];
my $outFile = $ARGV[3];

my %orthoPairs = ();
my %vogPairs = ();

open(FTW,">$outFile") or die;

open(FTR,"<$orthopairsFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    my $gNam1 = $arrInp[0];
    my $gNam2 = $arrInp[1];
    $orthoPairs{"$gNam1\t$gNam2"}=1;
}

close(FTR);

open(FTR,"<$vogsFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    my $gNam1 = $arrInp[0];
    my $gNam2 = $arrInp[1];
    $vogPairs{"$gNam1\t$gNam2"}=1;
    
}

close(FTR);


open(FTR,"<$identFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    my $gNam1 = $arrInp[0];
    my $gNam2 = $arrInp[1];
    
    if(not exists $vogPairs{"$gNam1\t$gNam2"})
    {
	if (not exists $orthoPairs{"$gNam1\t$gNam2"}) {
	    print FTW "$gNam1\t$gNam2\n";
	}
	
	
    }
    
}

close(FTR);




close(FTW);



