#!/usr/bin/perl

#my $type = "uds";
#my $type = "ccbar";
#my $type = "BpBm";
my $type = "B0B0bar";
#my $type = "pilnu";
#my $type = "pi0lnu";
#my $type = "etalnu";
#my $type = "etaplnu";
#my $type = "rhoClnu";
#my $type = "rho0lnu";
#my $type = "omegalnu";

#my $mode = "continuum";
my $mode = "genBB";
#my $mode = $type;

my $run = "Run4";
#my $run = "Run3";
#my $run = "Run2";
#my $run = "Run1";

my $NbMin = 601;
my $NbMax = 866;
my $nb = $NbMin;

my $BaseDir = "/nfs/slac/g/udem/u1/users/ProdV7";

print "generating file MakeNtuples-$type-$run-$NbMin-$NbMax.script \n";
open(OUT,">MakeNtuples-$type-$run-$NbMin-$NbMax.script") or die;
print OUT "\n";
print OUT "echo Producing ntuples for all the skimmed $type $run collections...\n";
print OUT "\n";
print OUT "setenv BetaModeString $mode\n"; 
print OUT "setenv Patch MC\n";
print OUT "setenv lookMC true\n";
print OUT "setenv BaseDir $BaseDir \n";
print OUT "unsetenv BetaColl\n";
print OUT "\n\n\n";

while($nb<=$NbMax) 
{
    print OUT "echo Sending job $type-$run-$nb... \n";	
    print OUT "setenv BetaCollFile $BaseDir/Skimmed_DATASETS/SP-$type-XSLBtoXulnuFilter-$run-R14-$nb.tcl \n";
    print OUT "setenv RootFile  $BaseDir/ntuples/$type-$run-$nb.root \n";
    print OUT "bsub -q long -o $BaseDir/logFiles/$type-$run-$nb.out ../bin/Linux24RH72_i386_gcc2953/BetaMiniApp ../XslReader/runReaderMC.tcl \n";
    print OUT "\n";
    ++$nb;
}

print OUT "echo Voila \n\n";
close OUT;






