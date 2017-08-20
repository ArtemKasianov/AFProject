use strict;
use OrthologGraphList;
use OrthologousGroups;

my $identFile = $ARGV[0];
my $predictExpressionFile = $ARGV[1];
my $predictCodonFreqFile = $ARGV[2];
my $trivialNamesFile = $ARGV[3];




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

$orthologGraph->InitalizePredictTable("predict.table");
$orthologGraph->RemoveEdgesWithWeightByPredictionTable();
$orthologGraph->PrintOrthogroups("orthogroups.list");


