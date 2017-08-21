use strict;

my $athExpressionFile = $ARGV[0];
my $fescExpressionFile = $ARGV[1];
my $orthologPairsFile = $ARGV[2];
my $typeOfSample = $ARGV[3];
my $isHeader = $ARGV[4];
my $outFile = $ARGV[5];
my $outFilePairs = $ARGV[6];


open(FTW,">$outFile") or die;
open(FTW_1,">$outFilePairs") or die;


my %athExpression = ();
my %fescExpression = ();
my %athExonLens = ();
my %fescExonLens = ();
my %athDomain = ();
my %fescDomain = ();
my %athShannon = ();
my %fescShannon = ();
my %athFescExonVectorDist = ();
my %athFescEuclidDist = ();
my %athFescIdent = ();


sub LoadDataFromTabularFile
{
    my $fileName = $_[0];
    my $ptrHash = $_[1];
    my $isHeader = $_[2];
    
    open(FTR,"<$fileName") or die("$fileName");
    
    if ($isHeader == 1) {
	<FTR>;
    }
    
    
    while (my $input = <FTR>) {
	chomp($input);
	
	my @arrInp = split(/\t/,$input);
	
	my $gNam = $arrInp[0];
	
	my @arrGNam = split(/\./,$gNam);
	$gNam = $arrGNam[0];
	
	my $expData = $arrInp[1];
	for(my $i = 2;$i <= $#arrInp;$i++)
	{
	    my $currExp = $arrInp[$i];
	    $expData = $expData."\t".$currExp;
	}
	$ptrHash->{"$gNam"} = $expData;
    }
    close(FTR);
    
    
}

sub LoadDataFromTabularFileWithPairs
{
    my $fileName = $_[0];
    my $ptrHash = $_[1];
    my $isHeader = $_[2];
    
    open(FTR,"<$fileName") or die;
    if ($isHeader == 1) {
	<FTR>;
    }
    while (my $input = <FTR>) {
	chomp($input);
	
	my @arrInp = split(/\t/,$input);
	
	my $gNam1 = $arrInp[0];
	my @arrGNam = split(/\./,$gNam1);
	$gNam1 = $arrGNam[0];
	my $gNam2 = $arrInp[1];
	@arrGNam = split(/\./,$gNam2);
	$gNam2 = $arrGNam[0];
	my $distVal = $arrInp[2];
	if (exists $ptrHash->{"$gNam1"}) {
	    $ptrHash->{"$gNam1"}->{"$gNam2"} = $distVal;
	}
	else
	{
	    my %hashTmp = ();
	    $hashTmp{"$gNam2"}=$distVal;
	    $ptrHash->{"$gNam1"} = \%hashTmp;
	}
	
    }
    close(FTR);
    
    
}

sub PrintValsFromString
{
    my $gNam = $_[0];
    my $ptrHash = $_[1];
    my $ptrIndex = $_[2];
    
    if (exists $ptrHash->{"$gNam"}) {
	
	my $strVal = $ptrHash->{"$gNam"};
	my @arrStr = split(/\t/,$strVal);
	
	for(my $i = 0;$i <= $#arrStr;$i++)
	{
	    my $currStrVal = $arrStr[$i];
	    $$ptrIndex = $$ptrIndex + 1;
	    print FTW "\t$$ptrIndex:$currStrVal";
	}
	#print FTW "\n";
    }
    else
    {
	die("$gNam\n");
    }
    
    
}


#sub PrintValsFromStringDomains
#{
#    my $gNam = $_[0];
#    my $ptrHash = $_[1];
#    my $ptrIndex = $_[2];
#    
#    if (exists $ptrHash->{"$gNam"}) {
#	
#	my $strVal = $ptrHash->{"$gNam"};
#	my @arrStr = split(/\t/,$strVal);
#	
#	for(my $i = 0;$i <= $#arrStr;$i++)
#	{
#	    my $currStrVal = $arrStr[$i];
#	    $$ptrIndex = $$ptrIndex + 1;
#	    print FTW "\t$$ptrIndex:$currStrVal";
#	}
#	#print FTW "\n";
#    }
#    else
#    {
#	for(my $i = 0;$i < $maxDomains;$i++ )
#	{
#	    print FTW "\t$$ptrIndex:0";
#	}
#    }
#    
#    
#}



sub PrintValsFromStringDomains
{
    my $gNam1 = $_[0];
    my $gNam2 = $_[1];
    my $ptrHash1 = $_[2];
    my $ptrHash2 = $_[3];
    my $ptrIndex = $_[4];
    my $gene1Number = 0;
    my $gene2Number = 0;
    my $distFlag = 0.5;
    if ((exists $ptrHash1->{"$gNam1"}) && (exists $ptrHash2->{"$gNam2"})) {
	my $strVal1 = $ptrHash1->{"$gNam1"};
	my @arrStr1 = split(/\t/,$strVal1);
	my $strVal2 = $ptrHash2->{"$gNam2"};
	my @arrStr2 = split(/\t/,$strVal2);
	$gene1Number = 0;
	$gene2Number = 0;
	for(my $i = 0;$i <= $#arrStr1;$i++)
	{
	    my $currStrVal1 = $arrStr1[$i];
	    my $currStrVal2 = $arrStr2[$i];
	    if ($currStrVal1 == 1) {
		$gene1Number++;
	    }
	    if ($currStrVal2 == 1) {
		$gene2Number++;
	    }
	    
	}
	
	my $dist2 = 0;
	for(my $i = 0;$i <= $#arrStr1;$i++)
	{
	    my $currStrVal1 = $arrStr1[$i];
	    my $currStrVal2 = $arrStr2[$i];
	    my $diff = ($currStrVal2 - $currStrVal1);
	    my $diff2 = $diff*$diff;
	    $dist2 = $dist2 + $diff2;
	}
	
	$distFlag = 1-($dist2/($gene1Number+$gene2Number));
	
	
    }
    
    if ((not exists $ptrHash1->{"$gNam1"}) && (exists $ptrHash2->{"$gNam2"})) {
	my $strVal2 = $ptrHash2->{"$gNam2"};
	my @arrStr2 = split(/\t/,$strVal2);
	$gene1Number = 0;
	$gene2Number = 0;
	for(my $i = 0;$i <= $#arrStr2;$i++)
	{
	    my $currStrVal2 = $arrStr2[$i];
	    if ($currStrVal2 == 1) {
		$gene2Number++;
	    }
	    
	}
	
	$distFlag = 0;
	
    }
    if ((exists $ptrHash1->{"$gNam1"}) && (not exists $ptrHash2->{"$gNam2"})) {
	my $strVal1 = $ptrHash1->{"$gNam1"};
	my @arrStr1 = split(/\t/,$strVal1);
	$gene1Number = 0;
	$gene2Number = 0;
	for(my $i = 0;$i <= $#arrStr1;$i++)
	{
	    my $currStrVal1 = $arrStr1[$i];
	    if ($currStrVal1 == 1) {
		$gene1Number++;
	    }
	    
	}
	$distFlag = 0;
    }
    
    if ((not exists $ptrHash1->{"$gNam1"}) && (not exists $ptrHash2->{"$gNam2"})) {
	$gene1Number = 0;
	$gene2Number = 0;
	$distFlag = 0.5;
    }
    
    $$ptrIndex = $$ptrIndex + 1;
    print FTW "\t$$ptrIndex:$gene1Number";
    $$ptrIndex = $$ptrIndex + 1;
    print FTW "\t$$ptrIndex:$gene2Number";
    $$ptrIndex = $$ptrIndex + 1;
    print FTW "\t$$ptrIndex:$distFlag";
}


sub PrintValsFromStringPairsIdentity
{
    my $gNam1 = $_[0];
    my $gNam2 = $_[1];
    my $ptrHash = $_[2];
    my $ptrIndex = $_[3];
    $$ptrIndex = $$ptrIndex + 1;
    if (exists $ptrHash->{"$gNam1"}->{"$gNam2"}) {
	
	
	my $distVal = $ptrHash->{"$gNam1"}->{"$gNam2"};
	print FTW "\t$$ptrIndex:$distVal";
    }
    else
    {
	print FTW "\t$$ptrIndex:0.0";
    }
    
    
}
sub PrintValsFromStringPairs
{
    my $gNam1 = $_[0];
    my $gNam2 = $_[1];
    my $ptrHash = $_[2];
    my $ptrIndex = $_[3];
    $$ptrIndex = $$ptrIndex + 1;
    if (exists $ptrHash->{"$gNam1"}->{"$gNam2"}) {
	
	my $distVal = $ptrHash->{"$gNam1"}->{"$gNam2"};
	print FTW "\t$$ptrIndex:$distVal";
    }
    else
    {
	die("$gNam1\t$gNam2\n");
    }
    
    
}

LoadDataFromTabularFile($athExpressionFile,\%athExpression,0);
LoadDataFromTabularFile($fescExpressionFile,\%fescExpression,0);
#LoadDataFromTabularFile($athExonLensFile,\%athExonLens,0);
#LoadDataFromTabularFile($fescExonLensFile,\%fescExonLens,0);
#LoadDataFromTabularFile($athDomainsFile,\%athDomain,0);
#LoadDataFromTabularFile($fescDomainsFile,\%fescDomain,0);

#LoadDataFromTabularFileWithPairs($athFescExonVectorDistFile,\%athFescExonVectorDist,0);
#LoadDataFromTabularFileWithPairs($athFescEuclidDistFile,\%athFescEuclidDist,0);
#LoadDataFromTabularFileWithPairs($athFescIdentFile,\%athFescIdent,0);

open(FTR, "<$orthologPairsFile") or die;
if ($isHeader == 1) {
    <FTR>;
}
while (my $input = <FTR>) {
    chomp($input);
    
    my @arrInp = split(/\t/,$input);
    my $athNam = $arrInp[0];
    my $fescNam = $arrInp[1];
    
    my $index = 0;
    
    
    next if(not exists $athExpression{"$athNam"});
    next if(not exists $fescExpression{"$fescNam"});
    #next if(not exists $athExonLens{"$athNam"});
    #next if(not exists $fescExonLens{"$fescNam"});
    #next if(not exists $athDomain{"$athNam"});
    #next if(not exists $fescDomain{"$fescNam"});
    #next if(not exists $athFescExonVectorDist{"$athNam"}->{"$fescNam"});
    #next if(not exists $athFescEuclidDist{"$athNam"}->{"$fescNam"});
    print FTW "$typeOfSample";
    
    PrintValsFromString($athNam,\%athExpression,\$index);
    PrintValsFromString($fescNam,\%fescExpression,\$index);
    #PrintValsFromString($athNam,\%athExonLens,\$index);
    #PrintValsFromString($fescNam,\%fescExonLens,\$index);
    #PrintValsFromStringPairsIdentity($athNam,$fescNam,\%athFescIdent,\$index);
    #PrintValsFromStringPairs($athNam,$fescNam,\%athFescExonVectorDist,\$index);
    #PrintValsFromStringPairs($athNam,$fescNam,\%athFescEuclidDist,\$index);
    #PrintValsFromStringDomains($athNam,$fescNam,\%athDomain,\%fescDomain,\$index);
    print FTW "\n";
    print FTW_1 "$athNam\t$fescNam\n";
}




close(FTR);
close(FTW);
close(FTW_1);
