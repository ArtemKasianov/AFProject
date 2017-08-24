use strict;
use OrthologGraphList;
use OrthologousGroups;

my $identFile = $ARGV[0];
my $predictExpressionFile = $ARGV[1];
my $predictCodonFreqFile = $ARGV[2];
my $trivialNamesFile = $ARGV[3];
my $outBeforeFile = $ARGV[4];
my $outAfterFile = $ARGV[5];
my $predictFile = $ARGV[6];





my %predictExpressionPairs = ();
my %predictCodonFreqPairs = ();
my %identPairs = ();

my %allAthHash = ();
my %allFescHash = ();


my %trivialNames = ();

my %distrOtrhogroupsBeforeSizes = ();
my %distrOtrhogroupsAfterSizes = ();



open(FTR,"<$trivialNamesFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    my $trivNam = $arrInp[1];
    
    $trivialNames{"$athNam"}=$trivNam;
    
}




close(FTR);



my $countOrthoPairsAll = 0;



open(FTR,"<$predictExpressionFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    my $fescNam = $arrInp[1];
    my $predictVal = sprintf("%.2f",$arrInp[2]);
    
    $allAthHash{"$athNam"} = 1;
    $allFescHash{"$fescNam"} = 1;
    
    
    $predictExpressionPairs{"$athNam\t$fescNam"}=$predictVal;
    
    
}


close(FTR);



open(FTR,"<$predictCodonFreqFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    my $fescNam = $arrInp[1];
    my $predictVal = sprintf("%.2f",$arrInp[2]);
    
    $allAthHash{"$athNam"} = 1;
    $allFescHash{"$fescNam"} = 1;
    
    
    $predictCodonFreqPairs{"$athNam\t$fescNam"}=$predictVal;
    
    
}


close(FTR);


my @allAth = keys %allAthHash;
my @allFesc = keys %allFescHash;

open(FTR,"<$identFile") or die;

while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    
    my $athNam = $arrInp[0];
    
    
    
    my $fescNam = $arrInp[1];
    
    my $identVal = sprintf("%.2f",$arrInp[2]);
    if (exists $predictExpressionPairs{"$athNam\t$fescNam"}) {
	$identPairs{"$athNam\t$fescNam"}=$identVal;
    }
    
    
}


close(FTR);



my $orthologGraph = OrthologGraphList->new();
for(my $i = 0;$i <= $#allAth;$i++)
{
    my $currAth = $allAth[$i];
    for(my $j = 0;$j <= $#allFesc;$j++)
    {
	my $currFesc = $allFesc[$j];
	
	
	next if(not exists $predictExpressionPairs{"$currAth\t$currFesc"});
	next if(not exists $predictCodonFreqPairs{"$currAth\t$currFesc"});
	next if(not exists $identPairs{"$currAth\t$currFesc"});
	my $predictExpressionVal = $predictExpressionPairs{"$currAth\t$currFesc"};
	my $predictCodonFreqVal = $predictCodonFreqPairs{"$currAth\t$currFesc"};
	my $identVal = $identPairs{"$currAth\t$currFesc"};
	
	$orthologGraph->AddEdge($currAth,$currFesc,$predictExpressionVal,$identVal,$predictCodonFreqVal);
	
	
    }
}

my %hashVisited = ();
my %hashEdgesVisited = ();
my $numOrtho = 0;
my @arrGenesSingletons = ();


for(my $i = 0;$i <= $#allAth;$i++ )
{
    my $currAthFirst = $allAth[$i];
    #print "$currAthFirst\n";
    next if(exists $hashVisited{"$currAthFirst"});
    
    $hashVisited{"$currAthFirst"} = 1;
    
    my %tempVisitedHash = ();
    my @arrVertexies = ();
    push @arrVertexies,$currAthFirst;
    $orthologGraph->GetListOfVertexiesInConnectedComponent($currAthFirst,\%hashVisited,\@arrVertexies);
    if ($#arrVertexies == 0) {
	my $singletonVertex = $arrVertexies[0];
	die if($singletonVertex ne $currAthFirst);
	push @arrGenesSingletons,$singletonVertex;
	$distrOtrhogroupsBeforeSizes{1} += 1;
	next;
    }

    if (exists $distrOtrhogroupsBeforeSizes{($#arrVertexies+1)}) {
	$distrOtrhogroupsBeforeSizes{($#arrVertexies+1)} = $distrOtrhogroupsBeforeSizes{($#arrVertexies+1)} + 1;
    }
    else
    {
	$distrOtrhogroupsBeforeSizes{($#arrVertexies+1)} = 1;
    }
    
    
    my $orthoGroupObj = OrthologousGroups->new();
    
    
    for(my $i = 0;$i <= $#arrVertexies;$i++)
    {
        my $gNam = $arrVertexies[$i];
        $orthoGroupObj->AddGene($gNam);
    }
    
    my $ptrAthArr = $orthoGroupObj->GetAthArr();
    my $ptrFescArr = $orthoGroupObj->GetFescArr();
    
    my $orthologGraphIn = OrthologGraphList->new();
    for(my $i = 0;$i <= $#$ptrAthArr;$i++)
    {
	my $currAth = $ptrAthArr->[$i];
	for(my $j = 0;$j <= $#$ptrFescArr;$j++)
	{
	    my $currFesc = $ptrFescArr->[$j];
	    
	    
	    next if(not exists $predictExpressionPairs{"$currAth\t$currFesc"});
	    next if(not exists $predictCodonFreqPairs{"$currAth\t$currFesc"});
	    next if(not exists $identPairs{"$currAth\t$currFesc"});
	    my $predictExpressionVal = $predictExpressionPairs{"$currAth\t$currFesc"};
	    my $predictCodonFreqVal = $predictCodonFreqPairs{"$currAth\t$currFesc"};
	    my $identVal = $identPairs{"$currAth\t$currFesc"};
	    
	    $orthologGraphIn->AddEdge($currAth,$currFesc,$predictExpressionVal,$identVal,$predictCodonFreqVal);
	    
	    
	}
    }
    $orthologGraphIn->InitalizePredictTable("$predictFile");
    $orthologGraphIn->RemoveEdgesWithWeightByPredictionTable();
    
    
    
    
    my %hashVisitedIn = ();
    my %hashEdgesVisitedIn = ();
    my @arrGenesSingletonsIn = ();
    my $indexCutOrthoGroup = 0;
    
    for(my $i = 0;$i <= $#$ptrAthArr;$i++)
    {
	my $currAth = $ptrAthArr->[$i];
	
	next if(exists $hashVisitedIn{"$currAth"});
	$hashVisitedIn{"$currAth"} = 1;
	%tempVisitedHash = ();
	@arrVertexies = ();
	
	push @arrVertexies,$currAth;
	$orthologGraphIn->GetListOfVertexiesInConnectedComponent($currAth,\%hashVisitedIn,\@arrVertexies);
	if ($#arrVertexies == 0) {
	    my $singletonVertex = $arrVertexies[0];
	    die if($singletonVertex ne $currAth);
	    push @arrGenesSingletonsIn,$singletonVertex;
	    $distrOtrhogroupsAfterSizes{1} += 1;
	    next;
	}
	
	if ($numOrtho == 883) {
	    if ($indexCutOrthoGroup == 2) {
		print "arrVertexies - $#arrVertexies\n";
	    }
	    
	}
	
	
	if (exists $distrOtrhogroupsAfterSizes{($#arrVertexies+1)}) {
	    $distrOtrhogroupsAfterSizes{($#arrVertexies+1)} = $distrOtrhogroupsAfterSizes{($#arrVertexies+1)} + 1;
	}
	else
	{
	    $distrOtrhogroupsAfterSizes{($#arrVertexies+1)} = 1;
	}
    }
    
    for(my $i = 0;$i <= $#$ptrFescArr;$i++)
    {
	my $currFesc = $ptrFescArr->[$i];
	
	next if(exists $hashVisitedIn{"$currFesc"});
	$hashVisitedIn{"$currFesc"} = 1;
	%tempVisitedHash = ();
	@arrVertexies = ();
	
	push @arrVertexies,$currFesc;
	$orthologGraphIn->GetListOfVertexiesInConnectedComponent($currFesc,\%hashVisitedIn,\@arrVertexies);
	if ($#arrVertexies == 0) {
	    my $singletonVertex = $arrVertexies[0];
	    die if($singletonVertex ne $currFesc);
	    push @arrGenesSingletonsIn,$singletonVertex;
	    $distrOtrhogroupsAfterSizes{1} += 1;
	    next;
	}
	
	die;
    }
    
	
    $numOrtho++;
    

}

for(my $i = 0;$i <= $#allFesc;$i++ )
{
    my $currFesc = $allFesc[$i];
    #print "$currFesc\n";
    next if(exists $hashVisited{"$currFesc"});
    
    $hashVisited{"$currFesc"} = 1;
    
    my %tempVisitedHash = ();
    my @arrVertexies = ();
    push @arrVertexies,$currFesc;
    $orthologGraph->GetListOfVertexiesInConnectedComponent($currFesc,\%hashVisited,\@arrVertexies);
    if ($#arrVertexies == 0) {
	my $singletonVertex = $arrVertexies[0];
	die if($singletonVertex ne $currFesc);
	push @arrGenesSingletons,$singletonVertex;
	$distrOtrhogroupsBeforeSizes{1} += 1;
	next;
    }
    die;
}

open(FTW_B,">$outBeforeFile") or die;

my @arrSizes = sort {$a <=> $b} keys %distrOtrhogroupsBeforeSizes;

for(my $i = 0;$i <= $#arrSizes;$i++)
{
    my $currSize = $arrSizes[$i];
    my $currCount = $distrOtrhogroupsBeforeSizes{$currSize};
    print FTW_B "$currSize\t$currCount\n";
}



close(FTW_B);


open(FTW_A,">$outAfterFile") or die;

@arrSizes = sort {$a <=> $b} keys %distrOtrhogroupsAfterSizes;

for(my $i = 0;$i <= $#arrSizes;$i++)
{
    my $currSize = $arrSizes[$i];
    my $currCount = $distrOtrhogroupsAfterSizes{$currSize};
    print FTW_A "$currSize\t$currCount\n";
}



close(FTW_A);



