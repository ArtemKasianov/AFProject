use strict;
use OrthologousGroups;


my $orthoFile = $ARGV[0];
my $outFile = $ARGV[1];


my %pairsUsed = ();


open(FTW,">$outFile") or die;

open(FTR,"<$orthoFile") or die;

my $lastPairs = "";

my $countOfPairs = 0;
while (my $input = <FTR>) {
    chomp($input);
    
    
    my $orthoGroups = OrthologousGroups->new();
    
    my @arrInp = split(/\t/,$input);
    
    my $orthoAthNam = $arrInp[0];
    my $orthoFescNam = $arrInp[1];
    
    for(my $i = 2;$i <= $#arrInp;$i++)
    {
	my $gNam = $arrInp[$i];
	$orthoGroups->AddGene($gNam);
	
    }
    
    
    my $ptrArrAth = $orthoGroups->GetAthArr();
    my $ptrArrFesc = $orthoGroups->GetFescArr();
    
    #print FTW "$orthoAthNam\t$orthoFescNam\n";
    
    for(my $i = 0;$i <= $#$ptrArrFesc;$i++)
    {
	my $currFesc = $ptrArrFesc->[$i];
	
	if (exists $pairsUsed{"$orthoAthNam\t$currFesc"}) {
	    next;
	}
	else
	{
	    print FTW "$orthoAthNam\t$currFesc\n";
	    $pairsUsed{"$orthoAthNam\t$currFesc"} = 1;
	}
	
	
	
    }
    
    for(my $i = 0;$i <= $#$ptrArrAth;$i++)
    {
	my $currAth = $ptrArrAth->[$i];
	if (exists $pairsUsed{"$currAth\t$orthoFescNam"}) {
	    next;
	}
	else
	{
	    print FTW "$currAth\t$orthoFescNam\n";
	    $pairsUsed{"$currAth\t$orthoFescNam"} = 1;
	}
    }
    
    
    
}



close(FTR);


close(FTW);
