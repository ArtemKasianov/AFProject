use strict;

my $athHitsFile = $ARGV[0];
my $fescHitsFile = $ARGV[1];
my $mRNAListAth = $ARGV[2];
my $orthoPairsFile = $ARGV[3];
my $outVOGFile = $ARGV[4];




my %athHits = ();
my %fescHits = ();
my %athmRNA = ();


open(FTW,">$outVOGFile") or die;


open(FTR,"<$mRNAListAth") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    $athmRNA{"$input"}=1;
    
}


close(FTR);


open(FTR,"<$athHitsFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    
    next if(not exists $athmRNA{"$athNam"});
    my @arrTmp = ();
    
    for(my $i = 1;$i <= $#arrInp;$i++)
    {
	my $currFesc = $arrInp[$i];
	
	push @arrTmp,$currFesc;
    }
    $athHits{"$athNam"} = \@arrTmp;
}

close(FTR);

open(FTR,"<$fescHitsFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    
    my $fescNam = $arrInp[0];
    my @arrTmp = ();
    
    for(my $i = 1;$i <= $#arrInp;$i++)
    {
	my $currAth = $arrInp[$i];
	next if(not exists $athmRNA{"$currAth"});
	push @arrTmp,$currAth;
    }
    $fescHits{"$fescNam"} = \@arrTmp;
}

close(FTR);

open(FTR,"<$orthoPairsFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    my @arrInp = split(/\t/,$input);
    my $athNam = $arrInp[0];
    my $fescNam = $arrInp[1];
    
    if(not exists $athHits{"$athNam"})
    {
	print "$athNam\n";
	next;
    }
    if(not exists $fescHits{"$fescNam"})
    {
	print "$fescNam\n";
	next;
    }
    
    print FTW "$athNam\t$fescNam";
    my $ptrArrAthHits = $athHits{"$athNam"};
    my $ptrArrFescHits = $fescHits{"$fescNam"};
    
    
    my $sizeOfVOG = $#$ptrArrAthHits + $#$ptrArrFescHits + 2;
    print "$sizeOfVOG\n";
    
    for(my $i = 0;$i <= $#$ptrArrAthHits;$i++)
    {
	my $currAth = $ptrArrAthHits->[$i];
	next if($currAth eq $fescNam);
	print FTW "\t$currAth";
    }
    
    for(my $i = 0;$i <= $#$ptrArrFescHits;$i++)
    {
	my $currFesc = $ptrArrFescHits->[$i];
	next if($currFesc eq $athNam);
	print FTW "\t$currFesc";
    }
    print FTW "\n";
    
}



close(FTR);


close(FTW);



