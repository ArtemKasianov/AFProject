use strict;




my $expressionFile = $ARGV[0];
my $athMrnaFile = $ARGV[1];
my $outFile = $ARGV[2];

my %mRNAList = ();


open(FTR,"<$athMrnaFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    $mRNAList{"$input"}=1;
}

close(FTR);


open(FTW,">$outFile") or die;


open(FTR,"<$expressionFile") or die;

my $title = <FTR>;
print FTW "$title";

while(my $input = <FTR>)
{
    chomp($input);
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    if (exists $mRNAList{"$athNam"}) {
        print FTW "$input\n";
    }
    
    
}

close(FTR);

close(FTW);