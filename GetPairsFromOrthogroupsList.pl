use strict;
use OrthologGraphList;
use OrthologousGroups;

my $orthogroupsListFile = $ARGV[0];
my $orthoMCLPairsFile = $ARGV[1];
my $outPairsFile = $ARGV[2];


my %pairsList = ();
my $countUnionList = 0;


open(FTR,"<$orthoMCLPairsFile") or die;


while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    my $fescNam = $arrInp[1];
    
    $pairsList{"$athNam\t$fescNam"} = 1;
    
}



close(FTR);

open(FTW,">$outPairsFile") or die;

open(FTR, "<$orthogroupsListFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    my $orthoGroupObj = OrthologousGroups->new();
    
    
    for(my $i = 0;$i <= $#arrInp;$i++)
    {
        my $gNam = $arrInp[$i];
        $orthoGroupObj->AddGene($gNam);
    }
    
    my $ptrAthArr = $orthoGroupObj->GetAthArr();
    my $ptrFescArr = $orthoGroupObj->GetFescArr();
    
    
    my $sizeAth = $#$ptrAthArr+1;
    my $sizeFesc = $#$ptrFescArr+1;
    
    my $sizeOfOrthoGroup = $sizeAth + $sizeFesc;
    next if($sizeOfOrthoGroup != 2);
    
    my $currAth = $ptrAthArr->[0];
    my $currFesc = $ptrFescArr->[0];
    print FTW "$currAth\t$currFesc\n";
    
    if (exists $pairsList{"$currAth\t$currFesc"}) {
	$countUnionList++;
    }
    
    
    
    
}

close(FTR);


close(FTW);


print "$countUnionList\n";

