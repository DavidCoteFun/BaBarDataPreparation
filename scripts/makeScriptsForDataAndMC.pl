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

my @runs = ("Run1","Run2","Run3","Run4");
foreach $run (@runs) {
    my @types = ("uds","ccbar","tautau","BpBm","B0B0bar","OnPeak","OffPeak");
    foreach $type (@types) {
	
	my $mode = "";
	if(($type eq "uds") || ($type eq "ccbar") || ($type eq "tautau")){ $mode="continuum"; }
	elsif(($type eq "BpBm") || ($type eq "B0B0bar")){ $mode="genBB"; }
	elsif(($type eq "OnPeak") || ($type eq "OffPeak")){ $mode="data"; }	

	
	if( !($mode eq "data") ) { 

	    $NbMax=0;
	    $NbMin=1;
	    my @files = `ls -1 ../Skimmed_DATASETS/*-$type-*-$run-*`;
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
		print OUT "setenv BetaModeString $mode\n"; 
		print OUT "setenv Patch MC\n";
		print OUT "setenv lookMC true\n";
		print OUT "setenv BaseDir $BaseDir \n";
		print OUT "unsetenv BetaColl\n";
		print OUT "\n\n\n";
		
		my $nb = $NbMin;
		while($nb<=$NbMaxTmp) 
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

		$NbMin+=$UltimateMax;
		$NbMax-=$UltimateMax;
		if($NbMax>$UltimateMax) { $NbMaxTmp+=$UltimateMax; }
		elsif(true){ $NbMaxTmp+=$NbMax; }
	    }	    
	}


	elsif($mode eq "data") {

	    $NbMax=0;
	    $NbMin=1;
	    my @files = `ls -1 ../Skimmed_DATASETS/*-$run-$type-*`;
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
		print OUT "setenv Patch Run2\n";
		print OUT "setenv BaseDir $BaseDir \n";
		print OUT "unsetenv BetaModeString\n"; 
		print OUT "unsetenv BetaColl\n";
		print OUT "unsetenv lookMC\n";
		print OUT "\n\n\n";
		
		my $nb = $NbMin;
		while($nb<=$NbMaxTmp) 
		{
		    print OUT "echo Sending job $type-$run-$nb... \n";
		    print OUT "setenv BetaCollFile $BaseDir/Skimmed_DATASETS/XSLBtoXulnuFilter-$run-$type-R14a-Total-$nb.tcl \n"; 
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
}














