use strict;
use OrthologGraphList;
use OrthologousGroups;

my $identFile = $ARGV[0];
my $predictExpressionFile = $ARGV[1];
my $predictCodonFreqFile = $ARGV[2];
my $trivialNamesFile = $ARGV[3];
my $outDir = $ARGV[4];
my $predictFile = $ARGV[5];




my %predictExpressionPairs = ();
my %predictCodonFreqPairs = ();
my %identPairs = ();

my %allAthHash = ();
my %allFescHash = ();


my %trivialNames = ();


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
    $orthologGraph->GetListOfVertexiesInConnectedComponent($currAthFirst,\%tempVisitedHash,\@arrVertexies);
    if ($#arrVertexies == 0) {
	my $singletonVertex = $arrVertexies[0];
	die if($singletonVertex ne $currAthFirst);
	push @arrGenesSingletons,$singletonVertex;
	next;
    }
    mkdir("$outDir/$numOrtho");
    open(FTW,">$outDir/$numOrtho/beforeCut.orthogroup.$numOrtho.graphviz") or die;
    print FTW "graph G {\n";
    close(FTW);
    $orthologGraph->PrintOrthogroupInFile("$currAthFirst",\%hashVisited,0,\%hashEdgesVisited,"$outDir/$numOrtho/beforeCut.orthogroup.$numOrtho.graphviz",\%trivialNames);
    open(FTW,">>$outDir/$numOrtho/beforeCut.orthogroup.$numOrtho.graphviz") or die;
    print FTW "}\n";
    
    close(FTW);
    
    
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
	$orthologGraphIn->GetListOfVertexiesInConnectedComponent($currAth,\%tempVisitedHash,\@arrVertexies);
	if ($#arrVertexies == 0) {
	    my $singletonVertex = $arrVertexies[0];
	    die if($singletonVertex ne $currAth);
	    push @arrGenesSingletonsIn,$singletonVertex;
	    next;
	}
	
	if ($numOrtho == 883) {
	    if ($indexCutOrthoGroup == 2) {
		print "arrVertexies - $#arrVertexies\t$currAth\n";
	    }
	    if ($indexCutOrthoGroup == 1) {
		print "arrVertexies - $#arrVertexies\t$currAth\n";
	    }
	    if ($indexCutOrthoGroup == 0) {
		print "arrVertexies - $#arrVertexies\t$currAth\n";
	    }
	    
	}
	
	
	open(FTW,">$outDir/$numOrtho/afterCut.orthogoroup.$numOrtho.$indexCutOrthoGroup.graphviz") or die;
	print FTW "graph G {\n";
	close(FTW);
	$orthologGraphIn->PrintOrthogroupInFile("$currAth",\%hashVisitedIn,0,\%hashEdgesVisitedIn,"$outDir/$numOrtho/afterCut.orthogoroup.$numOrtho.$indexCutOrthoGroup.graphviz",\%trivialNames);
	open(FTW,">>$outDir/$numOrtho/afterCut.orthogoroup.$numOrtho.$indexCutOrthoGroup.graphviz") or die;
	print FTW "}\n";
	$indexCutOrthoGroup++;
	close(FTW);
    }
    
    for(my $i = 0;$i <= $#$ptrFescArr;$i++)
    {
	my $currFesc = $ptrFescArr->[$i];
	
	next if(exists $hashVisitedIn{"$currFesc"});
	$hashVisitedIn{"$currFesc"} = 1;
	%tempVisitedHash = ();
	@arrVertexies = ();
	
	push @arrVertexies,$currFesc;
	$orthologGraphIn->GetListOfVertexiesInConnectedComponent($currFesc,\%tempVisitedHash,\@arrVertexies);
	if ($#arrVertexies == 0) {
	    my $singletonVertex = $arrVertexies[0];
	    die if($singletonVertex ne $currFesc);
	    push @arrGenesSingletonsIn,$singletonVertex;
	    next;
	}
	
	die;
    }
    if ($#arrGenesSingletonsIn > -1) {
	open(FTW,">$outDir/$numOrtho/orthogoroup.$numOrtho.singletons.txt") or die;
	
	for(my $i = 0;$i <= $#arrGenesSingletonsIn;$i++)
	{
	    my $currGeneSingleton = $arrGenesSingletonsIn[$i];
	    
	    
	    if (substr($currGeneSingleton,0,2) eq "AT") {
		my $trivNam = $currGeneSingleton;
		if (exists $trivialNames{"$currGeneSingleton"}) {
		    $trivNam = $trivialNames{"$currGeneSingleton"};
		}
		$currGeneSingleton = "$currGeneSingleton\_$trivNam";
	    }
	    
	    
	    print FTW "$currGeneSingleton\n";
	}
	
	close(FTW);
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
    $orthologGraph->GetListOfVertexiesInConnectedComponent($currFesc,\%tempVisitedHash,\@arrVertexies);
    if ($#arrVertexies == 0) {
	my $singletonVertex = $arrVertexies[0];
	die if($singletonVertex ne $currFesc);
	push @arrGenesSingletons,$singletonVertex;
	next;
    }
    die;
}

if ($#arrGenesSingletons > -1) {
    open(FTW,">$outDir/singletons.txt") or die;
    
    for(my $i = 0;$i <= $#arrGenesSingletons;$i++)
    {
	my $currGeneSingleton = $arrGenesSingletons[$i];
	
	
	if (substr($currGeneSingleton,0,2) eq "AT") {
	    my $trivNam = $currGeneSingleton;
	    if (exists $trivialNames{"$currGeneSingleton"}) {
		$trivNam = $trivialNames{"$currGeneSingleton"};
	    }
	    $currGeneSingleton = "$currGeneSingleton\_$trivNam";
	}
	
	
	print FTW "$currGeneSingleton\n";
    }
    
    close(FTW);
}

