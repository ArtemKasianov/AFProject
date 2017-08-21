use strict;




my $identFile = $ARGV[0];
my $treshold = $ARGV[1];
my $outFile = $ARGV[2];


open(FTW,">$outFile") or die;


open(FTR,"<$identFile") or die;

while(my $input = <FTR>)
{
    chomp($input);
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    my $fescNam = $arrInp[1];
    my $identVal = $arrInp[2];
    
    if($identFile >= $treshold)
    {
        print FTW "$input\n";
    }
    
}

close(FTR);

close(FTW);