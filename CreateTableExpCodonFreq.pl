use strict;

my $expressionFile = $ARGV[0];
my $codonFreqFile = $ARGV[1];
my $outFile = $ARGV[2];

my %expressionPairs = ();
my %codonFreqPairs = ();


open(FTW,">$outFile") or die;


open(FTR,"<$expressionFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    my $fescNam = $arrInp[1];
    my $expressionVal = $arrInp[2];
    
    
    $expressionPairs{"$athNam\t$fescNam"}=$expressionVal;
    
}

close(FTR);



open(FTR,"<$codonFreqFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    my $fescNam = $arrInp[1];
    my $codonFreqVal = $arrInp[2];
    
    
    $codonFreqPairs{"$athNam\t$fescNam"}=$codonFreqVal;
    
}

close(FTR);


my @arrPairs = keys %expressionPairs;

print FTW "AthNam\tFescNam\tExpression\tCodonFreq\n";
for(my $i = 0;$i <= $#arrPairs;$i++)
{
    my $currPair = $arrPairs[$i];
    
    die if(not exists $expressionPairs{"$currPair"});
    die if(not exists $codonFreqPairs{"$currPair"});
    
    my $expressionVal = $expressionPairs{"$currPair"};
    my $codonFreqVal = $codonFreqPairs{"$currPair"};
    print FTW "$currPair\t$expressionVal\t$codonFreqVal\n";
}



close(FTW);
