#!/usr/bin/perl

use Getopt::Std;
use strict;

my $type = "OnPeak";
#my $type = "OffPeak";

#my $run = "Run4";
#my $run = "Run3";
my $run = "Run2";
#my $run = "Run1";

my $NbMin = 368;
my $NbMax = 735;
my $nb = $NbMin;

my $BaseDir = "/nfs/slac/g/udem/u1/users/ProdV7";

print "generating file MakeNtuples-$type-$run-$NbMin-$NbMax.script \n";
open(OUT,">MakeNtuples-$type-$run-$NbMin-$NbMax.script") or die;
print OUT "\n";
print OUT "echo Producing ntuples for all the skimmed $type $run collections...\n";
print OUT "\n";
print OUT "setenv Patch Run2\n";
print OUT "setenv BaseDir $BaseDir \n";
print OUT "unsetenv BetaModeString\n"; 
print OUT "unsetenv BetaColl\n";
print OUT "unsetenv lookMC\n";
print OUT "\n\n\n";

while($nb<=$NbMax) 
{
    print OUT "echo Sending job $type-$run-$nb... \n";
    if(!($run eq "Run4")) { print OUT "setenv BetaCollFile $BaseDir/Skimmed_DATASETS/XSLBtoXulnuFilter-$run-$type-R14a-$nb.tcl \n"; }
    elsif(true){ print OUT "setenv BetaCollFile $BaseDir/Skimmed_DATASETS/XSLBtoXulnuFilter-$run-$type-R14-$nb.tcl \n"; }
    print OUT "setenv RootFile  $BaseDir/ntuples/$type-$run-$nb.root \n";
    print OUT "bsub -q long -o $BaseDir/logFiles/$type-$run-$nb.out ../bin/Linux24RH72_i386_gcc2953/BetaMiniApp ../XslReader/runReaderData.tcl \n";
    print OUT "\n";
    ++$nb;
}


print OUT "echo Voila \n\n";
close OUT;






