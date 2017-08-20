use strict;

my $tableFile = $ARGV[0];
my $outExpressionFile = $ARGV[1];
my $outIdentFile = $ARGV[2];
my $outCodonFreqFile = $ARGV[3];


open(FTW_E, ">$outExpressionFile") or die;
open(FTW_I, ">$outIdentFile") or die;
open(FTW_C, ">$outCodonFreqFile") or die;

open(FTR,"<$tableFile") or die;
<FTR>;
while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    my $fescNam = $arrInp[1];
    my $expVal = $arrInp[2];
    my $identVal = $arrInp[4];
    my $codonFreqVal = $arrInp[3];
    
    
    print FTW_E "$athNam\t$fescNam\t$expVal\n";
    print FTW_I "$athNam\t$fescNam\t$identVal\n";
    print FTW_C "$athNam\t$fescNam\t$codonFreqVal\n";
    
}





close(FTR);




close(FTW_E);
close(FTW_I);
close(FTW_C);


