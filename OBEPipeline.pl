use strict;

my $identTreshold = $ARGV[0];
my $identFile = $ARGV[1];
my $athExpressionFile = $ARGV[2];
my $fescExpressionFile = $ARGV[3];
my $athCodonFreqFile = $ARGV[4];
my $fescCodonFreqFile = $ARGV[5];
my $athMrnaFile = $ARGV[6];
my $orthopairsFile = $ARGV[7];
my $trivialNamesFile = $ARGV[8];


mkdir("results_treshold_$identTreshold");
mkdir("results_treshold_$identTreshold/data_for_learning");
system("perl FilterIdentFileByTreshold.pl $identFile $identTreshold results_treshold_$identTreshold/data_for_learning/Ath_vs_Fesc.gt$identTreshold.ident");
system("perl FilterIdentityByMRnaList.pl results_treshold_$identTreshold/data_for_learning/Ath_vs_Fesc.gt$identTreshold.ident $athMrnaFile results_treshold_$identTreshold/data_for_learning/Ath_vs_Fesc.gt$identTreshold.filtered.ident");
system("perl FilterExpressionFileByAthMrna.pl $athExpressionFile $athMrnaFile results_treshold_$identTreshold/data_for_learning/ath.expression.filtered.predict");
system("perl FilterCodonFreqFileByAthMrna.pl $athCodonFreqFile $athMrnaFile results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict");
system("perl FilterByMrnaAndIdentity.pl results_treshold_$identTreshold/data_for_learning/Ath_vs_Fesc.gt$identTreshold.filtered.ident $athMrnaFile $orthopairsFile results_treshold_$identTreshold/data_for_learning/orthopairs.filtered.list");
mkdir("results_treshold_$identTreshold/data_for_learning/folds");
mkdir("results_treshold_$identTreshold/data_for_learning/svm");
mkdir("results_treshold_$identTreshold/data_for_learning/virtual_OG");
system("perl GetRandomFoldOfOrthopairs.pl results_treshold_$identTreshold/data_for_learning/orthopairs.filtered.list 10 results_treshold_$identTreshold/data_for_learning/folds");
system("perl GetHitsOnlyForAthAndFesc.pl results_treshold_$identTreshold/data_for_learning/Ath_vs_Fesc.gt$identTreshold.filtered.ident results_treshold_$identTreshold/data_for_learning/virtual_OG/fesc_for_ath.list results_treshold_$identTreshold/data_for_learning/virtual_OG/ath_for_fesc.list temp.ath.stats temp.fesc.stats");

for(my $i = 0;$i <= 9;$i++)
{
    system("perl CreateVOGs.pl results_treshold_$identTreshold/data_for_learning/virtual_OG/fesc_for_ath.list results_treshold_$identTreshold/data_for_learning/virtual_OG/ath_for_fesc.list $athMrnaFile results_treshold_$identTreshold/data_for_learning/folds/fold_$i.orthopairs results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_$i.list");
    system("perl GetPairsFromOrthogroupsStr.pl results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_$i.list results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_$i.pairs.list");
}


mkdir("results_treshold_$identTreshold/training");
mkdir("results_treshold_$identTreshold/training/codonFreq");
for(my $i = 0;$i <= 9;$i++)
{
    mkdir("results_treshold_$identTreshold/training/codonFreq/folds_$i");
    mkdir("results_treshold_$identTreshold/training/codonFreq/folds_$i/model");
    
    my $strCat = "cat";
    
    for(my $j = 0;$j <= 9;$j++)
    {
	if ($i != $j) {
	    $strCat = $strCat." results_treshold_$identTreshold/data_for_learning/folds/fold_$j.orthopairs";
	}
    }
    $strCat = $strCat." >results_treshold_$identTreshold/training/codonFreq/folds_$i/not_fold_$i.orthopairs"; 
    system("$strCat");
    
    $strCat = "cat";
    
    for(my $j = 0;$j <= 9;$j++)
    {
	if ($i != $j) {
	    $strCat = $strCat." results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_$j.pairs.list";
	}
    }
    $strCat = $strCat." >results_treshold_$identTreshold/training/codonFreq/folds_$i/virtual.orthogroup.not_fold_$i.pairs.list"; 
    system("$strCat");
    system("perl GenerateSVMFile.codonFreq.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/training/codonFreq/folds_$i/not_fold_$i.orthopairs 1 0 results_treshold_$identTreshold/training/codonFreq/folds_$i/orthoMCL.oma.intersect.positive.examples.not_fold_$i.svm");
    system("perl GenerateSVMFile.codonFreq.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/training/codonFreq/folds_$i/virtual.orthogroup.not_fold_$i.pairs.list 0 0 results_treshold_$identTreshold/training/codonFreq/folds_$i/orthoMCL.oma.intersect.negative.examples.not_fold_$i.svm");
    system("perl GenerateSVMFile.codonFreq.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/data_for_learning/folds/fold_$i.orthopairs 1 0 results_treshold_$identTreshold/training/codonFreq/folds_$i/orthoMCL.oma.intersect.positive.examples.fold_$i.svm");
    system("perl GenerateSVMFile.codonFreq.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_$i.pairs.list 0 0 results_treshold_$identTreshold/training/codonFreq/folds_$i/orthoMCL.oma.intersect.negative.examples.fold_$i.svm");
    system("cat results_treshold_$identTreshold/training/codonFreq/folds_$i/orthoMCL.oma.intersect.positive.examples.not_fold_$i.svm results_treshold_$identTreshold/training/codonFreq/folds_$i/orthoMCL.oma.intersect.negative.examples.not_fold_$i.svm >results_treshold_$identTreshold/training/codonFreq/folds_$i/orthoMCL.oma.intersect.neg.pos.examples.not_fold_$i.svm");
    system("cat results_treshold_$identTreshold/training/codonFreq/folds_$i/orthoMCL.oma.intersect.negative.examples.fold_$i.svm results_treshold_$identTreshold/training/codonFreq/folds_$i/orthoMCL.oma.intersect.positive.examples.fold_$i.svm >results_treshold_$identTreshold/training/codonFreq/folds_$i/orthoMCL.oma.intersect.positive.negative.test.fold_$i.svm");
    system("python XGBoostTree.saveModel.py results_treshold_$identTreshold/training/codonFreq/folds_$i/orthoMCL.oma.intersect.neg.pos.examples.not_fold_$i.svm results_treshold_$identTreshold/training/codonFreq/folds_$i/orthoMCL.oma.intersect.positive.negative.test.fold_$i.svm 0 0.3 1 1 1 1 0 4 1 80 auc 20 results_treshold_$identTreshold/training/codonFreq/folds_$i/model/model_1_1_1_20.test");
    
}

mkdir("results_treshold_$identTreshold/training/codonFreq/folds_all");
mkdir("results_treshold_$identTreshold/training/codonFreq/folds_all/model");

my $strCat = "cat";

for(my $j = 0;$j <= 9;$j++)
{
    $strCat = $strCat." results_treshold_$identTreshold/data_for_learning/folds/fold_$j.orthopairs";
}
$strCat = $strCat." >results_treshold_$identTreshold/training/codonFreq/folds_all/fold_all.orthopairs"; 
system("$strCat");

$strCat = "cat";

for(my $j = 0;$j <= 9;$j++)
{
    $strCat = $strCat." results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_$j.pairs.list";
}
$strCat = $strCat." >results_treshold_$identTreshold/training/codonFreq/folds_all/virtual.orthogroup.fold_all.pairs.list"; 
system("$strCat");
system("perl GenerateSVMFile.codonFreq.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/training/codonFreq/folds_all/fold_all.orthopairs 1 0 results_treshold_$identTreshold/training/codonFreq/folds_all/orthoMCL.oma.intersect.positive.examples.fold_all.svm");
system("perl GenerateSVMFile.codonFreq.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/training/codonFreq/folds_all/virtual.orthogroup.fold_all.pairs.list 0 0 results_treshold_$identTreshold/training/codonFreq/folds_all/orthoMCL.oma.intersect.negative.examples.fold_all.svm");
system("perl GenerateSVMFile.codonFreq.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/data_for_learning/folds/fold_0.orthopairs 1 0 results_treshold_$identTreshold/training/codonFreq/folds_all/orthoMCL.oma.intersect.positive.examples.fold_0.svm");
system("perl GenerateSVMFile.codonFreq.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_0.pairs.list 0 0 results_treshold_$identTreshold/training/codonFreq/folds_all/orthoMCL.oma.intersect.negative.examples.fold_0.svm");
system("cat results_treshold_$identTreshold/training/codonFreq/folds_all/orthoMCL.oma.intersect.positive.examples.fold_all.svm results_treshold_$identTreshold/training/codonFreq/folds_all/orthoMCL.oma.intersect.negative.examples.fold_all.svm >results_treshold_$identTreshold/training/codonFreq/folds_all/orthoMCL.oma.intersect.neg.pos.examples.fold_all.svm");
system("cat results_treshold_$identTreshold/training/codonFreq/folds_all/orthoMCL.oma.intersect.negative.examples.fold_0.svm results_treshold_$identTreshold/training/codonFreq/folds_all/orthoMCL.oma.intersect.positive.examples.fold_0.svm >results_treshold_$identTreshold/training/codonFreq/folds_all/orthoMCL.oma.intersect.positive.negative.test.fold_0.svm");
system("python XGBoostTree.saveModel.py results_treshold_$identTreshold/training/codonFreq/folds_all/orthoMCL.oma.intersect.neg.pos.examples.fold_all.svm results_treshold_$identTreshold/training/codonFreq/folds_all/orthoMCL.oma.intersect.positive.negative.test.fold_0.svm 0 0.3 1 1 1 1 0 4 1 80 auc 20 results_treshold_$identTreshold/training/codonFreq/folds_all/model/model_1_1_1_20.test");






mkdir("results_treshold_$identTreshold/training");
mkdir("results_treshold_$identTreshold/training/expression");
for(my $i = 0;$i <= 9;$i++)
{
    mkdir("results_treshold_$identTreshold/training/expression/folds_$i");
    mkdir("results_treshold_$identTreshold/training/expression/folds_$i/model");
    
    my $strCat = "cat";
    
    for(my $j = 0;$j <= 9;$j++)
    {
	if ($i != $j) {
	    $strCat = $strCat." results_treshold_$identTreshold/data_for_learning/folds/fold_$j.orthopairs";
	}
    }
    $strCat = $strCat." >results_treshold_$identTreshold/training/expression/folds_$i/not_fold_$i.orthopairs"; 
    system("$strCat");
    
    $strCat = "cat";
    
    for(my $j = 0;$j <= 9;$j++)
    {
	if ($i != $j) {
	    $strCat = $strCat." results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_$j.pairs.list";
	}
    }
    $strCat = $strCat." >results_treshold_$identTreshold/training/expression/folds_$i/virtual.orthogroup.not_fold_$i.pairs.list"; 
    system("$strCat");
    system("perl GenerateSVMFile.expression.pl results_treshold_$identTreshold/data_for_learning/ath.expression.filtered.predict $fescExpressionFile results_treshold_$identTreshold/training/expression/folds_$i/not_fold_$i.orthopairs 1 0 results_treshold_$identTreshold/training/expression/folds_$i/orthoMCL.oma.intersect.positive.examples.not_fold_$i.svm");
    system("perl GenerateSVMFile.expression.pl results_treshold_$identTreshold/data_for_learning/ath.expression.filtered.predict $fescExpressionFile results_treshold_$identTreshold/training/expression/folds_$i/virtual.orthogroup.not_fold_$i.pairs.list 0 0 results_treshold_$identTreshold/training/expression/folds_$i/orthoMCL.oma.intersect.negative.examples.not_fold_$i.svm");
    system("perl GenerateSVMFile.expression.pl results_treshold_$identTreshold/data_for_learning/ath.expression.filtered.predict $fescExpressionFile results_treshold_$identTreshold/data_for_learning/folds/fold_$i.orthopairs 1 0 results_treshold_$identTreshold/training/expression/folds_$i/orthoMCL.oma.intersect.positive.examples.fold_$i.svm");
    system("perl GenerateSVMFile.expression.pl results_treshold_$identTreshold/data_for_learning/ath.expression.filtered.predict $fescExpressionFile results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_$i.pairs.list 0 0 results_treshold_$identTreshold/training/expression/folds_$i/orthoMCL.oma.intersect.negative.examples.fold_$i.svm");
    system("cat results_treshold_$identTreshold/training/expression/folds_$i/orthoMCL.oma.intersect.positive.examples.not_fold_$i.svm results_treshold_$identTreshold/training/expression/folds_$i/orthoMCL.oma.intersect.negative.examples.not_fold_$i.svm >results_treshold_$identTreshold/training/expression/folds_$i/orthoMCL.oma.intersect.neg.pos.examples.not_fold_$i.svm");
    system("cat results_treshold_$identTreshold/training/expression/folds_$i/orthoMCL.oma.intersect.negative.examples.fold_$i.svm results_treshold_$identTreshold/training/expression/folds_$i/orthoMCL.oma.intersect.positive.examples.fold_$i.svm >results_treshold_$identTreshold/training/expression/folds_$i/orthoMCL.oma.intersect.positive.negative.test.fold_$i.svm");
    system("python XGBoostTree.saveModel.py results_treshold_$identTreshold/training/expression/folds_$i/orthoMCL.oma.intersect.neg.pos.examples.not_fold_$i.svm results_treshold_$identTreshold/training/expression/folds_$i/orthoMCL.oma.intersect.positive.negative.test.fold_$i.svm 0 0.3 1 1 1 1 0 4 1 80 auc 20 results_treshold_$identTreshold/training/expression/folds_$i/model/model_1_1_1_20.test");
    
}

mkdir("results_treshold_$identTreshold/training/expression/folds_all");
mkdir("results_treshold_$identTreshold/training/expression/folds_all/model");

$strCat = "cat";

for(my $j = 0;$j <= 9;$j++)
{
    $strCat = $strCat." results_treshold_$identTreshold/data_for_learning/folds/fold_$j.orthopairs";
}
$strCat = $strCat." >results_treshold_$identTreshold/training/expression/folds_all/fold_all.orthopairs"; 
system("$strCat");

$strCat = "cat";

for(my $j = 0;$j <= 9;$j++)
{
    $strCat = $strCat." results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_$j.pairs.list";
}
$strCat = $strCat." >results_treshold_$identTreshold/training/expression/folds_all/virtual.orthogroup.fold_all.pairs.list"; 
system("$strCat");
system("perl GenerateSVMFile.codonFreq.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/training/expression/folds_all/fold_all.orthopairs 1 0 results_treshold_$identTreshold/training/expression/folds_all/orthoMCL.oma.intersect.positive.examples.fold_all.svm");
system("perl GenerateSVMFile.codonFreq.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/training/expression/folds_all/virtual.orthogroup.fold_all.pairs.list 0 0 results_treshold_$identTreshold/training/expression/folds_all/orthoMCL.oma.intersect.negative.examples.fold_all.svm");
system("perl GenerateSVMFile.codonFreq.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/data_for_learning/folds/fold_0.orthopairs 1 0 results_treshold_$identTreshold/training/expression/folds_all/orthoMCL.oma.intersect.positive.examples.fold_0.svm");
system("perl GenerateSVMFile.codonFreq.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_0.pairs.list 0 0 results_treshold_$identTreshold/training/expression/folds_all/orthoMCL.oma.intersect.negative.examples.fold_0.svm");
system("cat results_treshold_$identTreshold/training/expression/folds_all/orthoMCL.oma.intersect.positive.examples.fold_all.svm results_treshold_$identTreshold/training/expression/folds_all/orthoMCL.oma.intersect.negative.examples.fold_all.svm >results_treshold_$identTreshold/training/expression/folds_all/orthoMCL.oma.intersect.neg.pos.examples.fold_all.svm");
system("cat results_treshold_$identTreshold/training/expression/folds_all/orthoMCL.oma.intersect.negative.examples.fold_0.svm results_treshold_$identTreshold/training/expression/folds_all/orthoMCL.oma.intersect.positive.examples.fold_0.svm >results_treshold_$identTreshold/training/expression/folds_all/orthoMCL.oma.intersect.positive.negative.test.fold_0.svm");
system("python XGBoostTree.saveModel.py results_treshold_$identTreshold/training/expression/folds_all/orthoMCL.oma.intersect.neg.pos.examples.fold_all.svm results_treshold_$identTreshold/training/expression/folds_all/orthoMCL.oma.intersect.positive.negative.test.fold_0.svm 0 0.3 1 1 1 1 0 4 1 80 auc 20 results_treshold_$identTreshold/training/expression/folds_all/model/model_1_1_1_20.test");

$strCat = "cat";

for(my $j = 0;$j <= 9;$j++)
{
    $strCat = $strCat." results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_$j.pairs.list";
}
$strCat = $strCat." >results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_all.pairs.list";
system("$strCat");

system("perl GetPairsNotInOrthopairsAndVOGs.pl results_treshold_$identTreshold/data_for_learning/Ath_vs_Fesc.gt$identTreshold.filtered.ident results_treshold_$identTreshold/data_for_learning/orthopairs.filtered.list results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_all.pairs.list results_treshold_$identTreshold/data_for_learning/others.pairs");


mkdir("results_treshold_$identTreshold/predictions");
for(my $j = 0;$j <= 9;$j++)
{
    system("perl GenerateSVMFile.expression.pairs.pl results_treshold_$identTreshold/data_for_learning/ath.expression.filtered.predict $fescExpressionFile results_treshold_$identTreshold/data_for_learning/folds/fold_$j.orthopairs 0 0 results_treshold_$identTreshold/data_for_learning/svm/fold_$j.expression.svm results_treshold_$identTreshold/data_for_learning/svm/fold_$j.expression.pairs");
    system("perl GenerateSVMFile.expression.pairs.pl results_treshold_$identTreshold/data_for_learning/ath.expression.filtered.predict $fescExpressionFile results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_$j.pairs.list 0 0 results_treshold_$identTreshold/data_for_learning/svm/virtual.orthogroup.fold_$j.expression.svm results_treshold_$identTreshold/data_for_learning/svm/virtual.orthogroup.fold_$j.expression.pairs");
    system("python PredictByModelXGBoost.10.py results_treshold_$identTreshold/data_for_learning/svm/fold_$j.expression.svm results_treshold_$identTreshold/data_for_learning/svm/fold_$j.expression.pairs results_treshold_$identTreshold/training/expression/fold_$j/model/model_1_1_1_20.test 80 >results_treshold_$identTreshold/predictions/fold_$j.positive.orthopairs.expression.predictions");
    system("python PredictByModelXGBoost.10.py results_treshold_$identTreshold/data_for_learning/svm/virtual.orthogroup.fold_$j.expression.svm results_treshold_$identTreshold/data_for_learning/svm/virtual.orthogroup.fold_$j.expression.pairs results_treshold_$identTreshold/training/expression/fold_$j/model/model_1_1_1_20.test 80 >results_treshold_$identTreshold/predictions/fold_$j.negative.orthopairs.expression.predictions");
    
}

system("perl GenerateSVMFile.expression.pairs.pl results_treshold_$identTreshold/data_for_learning/ath.expression.filtered.predict $fescExpressionFile others.pairs 0 0 results_treshold_$identTreshold/data_for_learning/svm/others.expression.svm results_treshold_$identTreshold/data_for_learning/svm/others.expression.pairs");
system("python PredictByModelXGBoost.10.py results_treshold_$identTreshold/data_for_learning/svm/others.expression.svm results_treshold_$identTreshold/data_for_learning/svm/others.expression.pairs results_treshold_$identTreshold/training/expression/folds_all/model/model_1_1_1_20.test 80 >results_treshold_$identTreshold/predictions/others.expression.predictions");




for(my $j = 0;$j <= 9;$j++)
{
    system("perl GenerateSVMFile.expression.pairs.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/data_for_learning/folds/fold_$j.orthopairs 0 0 results_treshold_$identTreshold/data_for_learning/svm/fold_$j.codonFreq.svm results_treshold_$identTreshold/data_for_learning/svm/fold_$j.codonFreq.pairs");
    system("perl GenerateSVMFile.expression.pairs.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile results_treshold_$identTreshold/data_for_learning/virtual_OG/virtual.orthogroup.fold_$j.pairs.list 0 0 results_treshold_$identTreshold/data_for_learning/svm/virtual.orthogroup.fold_$j.codonFreq.svm results_treshold_$identTreshold/data_for_learning/svm/virtual.orthogroup.fold_$j.codonFreq.pairs");
    system("python PredictByModelXGBoost.10.py results_treshold_$identTreshold/data_for_learning/svm/fold_$j.codonFreq.svm results_treshold_$identTreshold/data_for_learning/svm/fold_$j.codonFreq.pairs results_treshold_$identTreshold/training/codonFreq/fold_$j/model/model_1_1_1_20.test 80 >results_treshold_$identTreshold/predictions/fold_$j.positive.orthopairs.codonFreq.predictions");
    system("python PredictByModelXGBoost.10.py results_treshold_$identTreshold/data_for_learning/svm/virtual.orthogroup.fold_$j.codonFreq.svm results_treshold_$identTreshold/data_for_learning/svm/virtual.orthogroup.fold_$j.codonFreq.pairs results_treshold_$identTreshold/training/codonFreq/fold_$j/model/model_1_1_1_20.test 80 >results_treshold_$identTreshold/predictions/fold_$j.negative.orthopairs.codonFreq.predictions");
    
}

system("perl GenerateSVMFile.expression.pairs.pl results_treshold_$identTreshold/data_for_learning/ath.codonFreq.filtered.predict $fescCodonFreqFile others.pairs 0 0 results_treshold_$identTreshold/data_for_learning/svm/others.codonFreq.svm results_treshold_$identTreshold/data_for_learning/svm/others.codonFreq.pairs");
system("python PredictByModelXGBoost.10.py results_treshold_$identTreshold/data_for_learning/svm/others.codonFreq.svm results_treshold_$identTreshold/data_for_learning/svm/others.codonFreq.pairs results_treshold_$identTreshold/training/codonFreq/folds_all/model/model_1_1_1_20.test 80 >results_treshold_$identTreshold/predictions/others.codonFreq.predictions");

$strCat = "cat";
my $strCatPositive = "cat";
my $strCatNegative = "cat";
for(my $j = 0;$j <= 9;$j++)
{
    $strCat = $strCat." results_treshold_$identTreshold/predictions/fold_$j.positive.orthopairs.codonFreq.predictions";
    $strCat = $strCat." results_treshold_$identTreshold/predictions/fold_$j.negative.orthopairs.codonFreq.predictions";
    
    $strCatPositive = $strCatPositive." results_treshold_$identTreshold/predictions/fold_$j.positive.orthopairs.codonFreq.predictions";
    $strCatNegative = $strCatNegative." results_treshold_$identTreshold/predictions/fold_$j.negative.orthopairs.codonFreq.predictions";
}
$strCat = $strCat." results_treshold_$identTreshold/predictions/others.codonFreq.predictions >results_treshold_$identTreshold/predictions/codonFreq.predictions";
$strCatPositive = $strCatPositive." >results_treshold_$identTreshold/predictions/positive.codonFreq.predictions";
$strCatNegative = $strCatNegative." >results_treshold_$identTreshold/predictions/negative.codonFreq.predictions";
system("$strCat");
system("$strCatPositive");
system("$strCatNegative");

$strCat = "cat";
$strCatPositive = "cat";
$strCatNegative = "cat";
for(my $j = 0;$j <= 9;$j++)
{
    $strCatPositive = $strCatPositive." results_treshold_$identTreshold/predictions/fold_$j.positive.orthopairs.expression.predictions";
    $strCatNegative = $strCatNegative." results_treshold_$identTreshold/predictions/fold_$j.negative.orthopairs.expression.predictions";
}
$strCat = $strCat." results_treshold_$identTreshold/predictions/others.expression.predictions >results_treshold_$identTreshold/predictions/expression.predictions";
$strCatPositive = $strCatPositive." >results_treshold_$identTreshold/predictions/positive.expression.predictions";
$strCatNegative = $strCatNegative." >results_treshold_$identTreshold/predictions/negative.expression.predictions";

system("$strCat");
system("$strCatPositive");
system("$strCatNegative");

system("perl CreateTableExpCodonFreq.pl results_treshold_$identTreshold/predictions/positive.expression.predictions results_treshold_$identTreshold/predictions/positive.codonFreq.predictions results_treshold_$identTreshold/predictions/positive.expressionCodonFreq.predictions");
system("perl CreateTableExpCodonFreq.pl results_treshold_$identTreshold/predictions/negative.expression.predictions results_treshold_$identTreshold/predictions/negative.codonFreq.predictions results_treshold_$identTreshold/predictions/negative.expressionCodonFreq.predictions");
system("python svm-pipe-part1.py results_treshold_$identTreshold/predictions/positive.expressionCodonFreq.predictions results_treshold_$identTreshold/predictions/negative.expressionCodonFreq.predictions results_treshold_$identTreshold/predictions/treshold.table");

mkdir("results_treshold_$identTreshold/graphs");
mkdir("results_treshold_$identTreshold/graphs/graphviz");
system("perl PrintOrthogroupsAfterCutByTreshold.beforeAfterCut.pl results_treshold_$identTreshold/data_for_learning/Ath_vs_Fesc.gt$identTreshold.ident results_treshold_$identTreshold/predictions/expression.predictions results_treshold_$identTreshold/predictions/codonFreq.predictions $trivialNamesFile results_treshold_$identTreshold/graphs/graphviz");
system("perl CountOrthogroupsDistrAfterCutByTreshold.beforeAfterCut.pl results_treshold_$identTreshold/data_for_learning/Ath_vs_Fesc.gt$identTreshold.ident results_treshold_$identTreshold/predictions/expression.predictions results_treshold_$identTreshold/predictions/codonFreq.predictions $trivialNamesFile results_treshold_$identTreshold/graphs/beforeCut.stats results_treshold_$identTreshold/graphs/afterCut.stats");

