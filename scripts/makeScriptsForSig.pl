#!/usr/bin/perl
use Getopt::Std;                                    # part of perl
use strict;

my $NbMin = 1;
my $NbMax = 0;
my $UltimateMax = 200;

my $BaseDir = "/nfs/slac/g/udem/u1/users/ProdV8";
my $run="";
my $type="";
my $file="";

my @runs = ("SP5","SP6");
foreach $run (@runs) {
    my @types = ("pilnu","pi0lnu","etalnu","etaplnu","rhoClnu","rho0lnu","omegalnu");
    foreach $type (@types) {
	
	my $SPnum = 0;
	if($type eq "pilnu") { $SPnum=4764; }
	elsif($type eq "pi0lnu") { $SPnum=4762; }
	elsif($type eq "etalnu") { $SPnum=4759; }
	elsif($type eq "etaplnu") { $SPnum=4760; }
	elsif($type eq "rhoClnu") { $SPnum=4766; }
	elsif($type eq "rho0lnu") { $SPnum=4763; }
	elsif($type eq "omegalnu") { $SPnum=4761; }
	

	$NbMax=0;
	$NbMin=1;
	my @files = `ls -1 ../Skimmed_DATASETS/SP-$SPnum-XSLBtoXulnuFilter-$run-*`;
	foreach $file (@files) {
	    ++$NbMax;
	}
	--$NbMax;
	
	my $NbMaxTmp=0;
	if($NbMax>$UltimateMax) { $NbMaxTmp=$UltimateMax; }
	elsif(true){ $NbMaxTmp=$NbMax; }
	while($NbMax>0) {
	    
	    print "generating file MakeNtuples-$type-$run-$NbMin-$NbMaxTmp.script \n";
	    open(OUT,">MakeNtuples-$type-$run-$NbMin-$NbMaxTmp.script") or die;
	    print OUT "\n";
	    print OUT "echo Producing ntuples for all the skimmed $type $run collections...\n";
	    print OUT "\n";
	    print OUT "setenv Patch MC\n";
	    print OUT "setenv lookMC true\n";
	    print OUT "setenv BaseDir $BaseDir \n";
	    print OUT "unsetenv BetaColl\n";
	    print OUT "setenv BetaModeString $type \n"; 
	    print OUT "\n\n\n";
	    
	    my $nb = $NbMin;
	    while($nb<=$NbMaxTmp) 
	    {
		print OUT "echo Sending job $type-$run-$nb... \n";
		print OUT "setenv BetaCollFile $BaseDir/Skimmed_DATASETS/SP-$SPnum-XSLBtoXulnuFilter-$run-$nb.tcl \n"; 
		print OUT "setenv RootFile  $BaseDir/ntuples/$type-$run-$nb.root \n";
		print OUT "bsub -q long -o $BaseDir/logFiles/$type-$run-$nb.out ../bin/Linux24RH72_i386_gcc2953/BetaMiniApp ../XslReader/runReaderData.tcl \n";
		print OUT "\n";
		++$nb;
	    }
	    
	    print OUT "echo Voila \n\n";
	    close OUT;
	    
	    $NbMin+=$UltimateMax;
	    $NbMax-=$UltimateMax;
	    if($NbMax>$UltimateMax) { $NbMaxTmp+=$UltimateMax; }
	    elsif(true){ $NbMaxTmp+=$NbMax; }		
	}	
    }		    
}














