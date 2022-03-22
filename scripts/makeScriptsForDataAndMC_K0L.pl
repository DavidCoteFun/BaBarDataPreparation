#!/usr/bin/perl
use Getopt::Std;                                    # part of perl
use strict;

my $UltimateMax = 200;
my $BaseDir = "/nfs/farm/babar/AWG26/cote/ProdV1";

my $run="";
my $type="";
my $file="";

my @runs = ("Run1","Run2","Run3","Run4");
foreach $run (@runs) {
    my @types = ("SP-5735","SP-5736","SP-1135","SP-1136","SP-1005","OnPeak","OffPeak");
    foreach $type (@types) {
	
	my $dumpMC = "true";
	if(($type eq "OnPeak") || ($type eq "OffPeak")){ $dumpMC="false"; }	

	#Get file name and nb of files
	my $NbMax=0;
	my $NbMin=1;
	my $name="";
	if($type eq "SP-1005") { 
	    $name="$type-AllEventsSkim-$run-R16a";
	    my @files = `ls -1 ../Collections/$name-*`;
	    foreach $file (@files) {
		++$NbMax;
	    }
	}	    
	elsif(($type eq "SP-1135") || ($type eq "SP-1136") ||($type eq "SP-5735") ||($type eq "SP-5736") ) { 
	    $name="$type-$run";
	    my @files = `ls -1 ../Collections/$name-*`;
	    foreach $file (@files) {
		++$NbMax;
	    }
	}	
	elsif(($type eq "OffPeak") || ($type eq "OnPeak") ) { 
	    $name="AllEventsSkim-$run-$type-R16a";
	    my @files = `ls -1 ../Collections/$name-*`;
	    foreach $file (@files) {
		++$NbMax;
	    }
	}	
	
	#write the scripts (to be sourced manually)
	my $NbMaxTmp=0;
	if($NbMax>$UltimateMax) { $NbMaxTmp=$UltimateMax; }
	elsif(true){ $NbMaxTmp=$NbMax; }
	while($NbMax>0) {

	    print "generating file MakeNtuples-$type-$run-$NbMin-$NbMaxTmp.script \n";
	    open(OUT,">MakeNtuples-$type-$run-$NbMin-$NbMaxTmp.script") or die;
	    print OUT "\n";
	    print OUT "echo Producing ntuples for all $type $run collections...\n";
	    print OUT "\n";
	    print OUT "unsetenv BetaColl\n";
	    print OUT "setenv DumpMCInfo $dumpMC\n";
	    print OUT "\n\n\n";
	    
	    my $nb = $NbMin;
	    while($nb<=$NbMaxTmp) 
	    {
		print OUT "echo Sending job $type-$run-$nb... \n";	
		print OUT "setenv BetaCollFile $BaseDir/Collections/$name-$nb.tcl \n";
		print OUT "setenv RootFile  $BaseDir/ntuples/$type-$run-$nb.root \n";
		print OUT "bsub -q xlong -o $BaseDir/logFiles/$type-$run-$nb.out ../bin/Linux24SL3_i386_gcc323-noOptimize-Debug/K0LMiniUser ../K0LUser/MyModule.tcl \n";
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














