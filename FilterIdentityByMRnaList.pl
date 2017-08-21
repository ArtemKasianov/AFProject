use strict;




my $identFile = $ARGV[0];
my $athMrnaFile = $ARGV[1];
my $outFile = $ARGV[2];

my %athMrnaList = ();

open(FTW,">$outFile") or die;

open(FTR,"<$athMrnaFile") or die;

while(my $input = <FTR>)
{
    chomp($input);
    
    $athMrnaList{"$input"} = 1;
    
}

close(FTR);


open(FTR,"<$identFile") or die;

while(my $input = <FTR>)
{
    chomp($input);
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    next if(not exists $athMrnaList{"$athNam"});
    
    print FTW "$input\n";
    
}

close(FTR);



close(FTW);
