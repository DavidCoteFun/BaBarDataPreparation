#!/usr/bin/perl


my $NbMin = 1;
my $NbMax = 0;

my $BaseDir = "/nfs/slac/g/udem/u1/users/ProdV7";


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
	    my @files = `ls -1 ../Skimmed_DATASETS/*-$type-*-$run-*`;
	    foreach $file (@files) {
		++$NbMax;
	    }
	    --$NbMax;


	    $myNbBad=0;
	
	    my $nb = $NbMin;
	    while($nb<=$NbMax) 
	    {
		my @fff = `grep -L fully $type-$run-$nb.out`;	
		foreach $fff (@fff) { print "$fff"; }
		++$nb;
	    }
	    print OUT "echo Voila \n\n";
	    close OUT;
	}


	elsif($mode eq "data") {

	    $NbMax=0;
	    my @files = `ls -1 ../Skimmed_DATASETS/*-$run-$type-*`;
	    foreach $file (@files) {
		++$NbMax;
	    }
	    --$NbMax;
	    
	    $myNbBad=0;
	    my $nb = $NbMin;
	    while($nb<=$NbMax) 
	    {
		my @fff = `grep -L fully $type-$run-$nb.out`;	
		foreach $fff (@fff) { print "$fff"; }
		++$nb;
	    }
	    
	}
	
	
    }
}














