use strict;

my $identFile = $ARGV[0];
my $outFileAth = $ARGV[1];
my $outFileFesc = $ARGV[2];
my $outFileAthStats = $ARGV[3];
my $outFileFescStats = $ARGV[4];



my %athHits = ();
my %fescHits = ();


open(FTW_ATH,">$outFileAth") or die;
open(FTW_FESC,">$outFileFesc") or die;
open(FTW_ATH_STATS,">$outFileAthStats") or die;
open(FTW_FESC_STATS,">$outFileFescStats") or die;



open(FTR,"<$identFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    my $fescNam = $arrInp[1];
    my $ident = $arrInp[2];
    
    
    if (exists $athHits{"$athNam"}) {
	my $ptrArrTmp = $athHits{"$athNam"};
	push @$ptrArrTmp,$fescNam;
    }
    else
    {
	my @arrTmp = ();
	push @arrTmp,$fescNam;
	$athHits{"$athNam"} = \@arrTmp;
    }
    
    if (exists $fescHits{"$fescNam"}) {
	my $ptrArrTmp = $fescHits{"$fescNam"};
	push @$ptrArrTmp,$athNam;
    }
    else
    {
	my @arrTmp = ();
	push @arrTmp,$athNam;
	$fescHits{"$fescNam"} = \@arrTmp;
    }
    
}





close(FTR);



my @arrAth = keys %athHits;

for(my $i = 0;$i <= $#arrAth;$i++)
{
    my $currAth = $arrAth[$i];
    my $ptrArrFescTmp = $athHits{"$currAth"};
    print FTW_ATH "$currAth";
    my $numberOfHits = $#$ptrArrFescTmp+1;
    print FTW_ATH_STATS "$currAth\t$numberOfHits\n";
    for(my $j = 0;$j <= $#$ptrArrFescTmp;$j++)
    {
	my $currFesc = $ptrArrFescTmp->[$j];
	print FTW_ATH "\t$currFesc";
    }
    print FTW_ATH "\n";
    
}



close(FTW_ATH);
close(FTW_ATH_STATS);




my @arrFesc = keys %fescHits;

for(my $i = 0;$i <= $#arrFesc;$i++)
{
    my $currFesc = $arrFesc[$i];
    my $ptrArrAthTmp = $fescHits{"$currFesc"};
    print FTW_FESC "$currFesc";
    my $numberOfHits = $#$ptrArrAthTmp+1;
    print FTW_FESC_STATS "$currFesc\t$numberOfHits\n";
    for(my $j = 0;$j <= $#$ptrArrAthTmp;$j++)
    {
	my $currAth = $ptrArrAthTmp->[$j];
	print FTW_FESC "\t$currAth";
    }
    print FTW_FESC "\n";
    
}



close(FTW_FESC_STATS);
close(FTW_FESC);




