use strict;


my $identPairsFile = $ARGV[0];
my $mRNAFile = $ARGV[1];
my $orthoPairsFile = $ARGV[2];
my $outFile = $ARGV[3];

my %identPairs = ();
my %mRNA = ();

open(FTW,">$outFile") or die;

open(FTR,"<$identPairsFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    
    my $fescNam = $arrInp[1];
    
    $identPairs{"$athNam\t$fescNam"} = 1;
    
    
}

close(FTR);


open(FTR,"<$mRNAFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    $mRNA{"$input"}=1;
    
}

close(FTR);



open(FTR,"<$orthoPairsFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    my @arrInp = split(/\t/,$input);
    my $athNam = $arrInp[0];
    my $fescNam = $arrInp[1];
    
    next if(not exists $identPairs{"$athNam\t$fescNam"});
    next if(not exists $mRNA{"$athNam"});
    
    print FTW "$input\n";
    
}

close(FTR);

close(FTW);


