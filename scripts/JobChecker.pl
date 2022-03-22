#!/usr/local/bin/perl

#scope of this script:
#        -grep log files to find out the ones that failed
#        -change the name of these file.out to file.out-Error
#        -write scripts (to be sourced manually) to rerun these jobs

my $BaseDir = "/nfs/slac/g/udem/u1/users/ProdV7";

#First create the head of the script files:
#continuum
print "generating file ../scripts/RerunJobs-Continuum.script \n";
open(cont,">../scripts/RerunJobs-Continuum.script") or die;
print cont "\n";
print cont "echo Producing ntuples for all the skimmed Continuum collections that once failed...\n";
print cont "\n";
print cont "setenv BetaModeString continuum\n"; 
print cont "setenv Patch MC\n";
print cont "setenv lookMC true\n";
print cont "setenv BaseDir $BaseDir \n";
print cont "unsetenv BetaColl\n";
print cont "\n\n\n";

#genBB
print "generating file ../scripts/RerunJobs-genBB.script \n";
open(BB,">../scripts/RerunJobs-genBB.script") or die;
print BB "\n";
print BB "echo Producing ntuples for all the skimmed genBB collections that once failed...\n";
print BB "\n";
print BB "setenv BetaModeString genBB\n"; 
print BB "setenv Patch MC\n";
print BB "setenv lookMC true\n";
print BB "setenv BaseDir $BaseDir \n";
print BB "unsetenv BetaColl\n";
print BB "\n\n\n";

#data
print "generating file ../scripts/RerunJobs-data.script \n";
open(data,">../scripts/RerunJobs-data.script") or die;
print data "\n";
print data "echo Producing ntuples for all the skimmed data collections that once failed...\n";
print data "\n";
print data "setenv Patch Run2\n";
print data "setenv BaseDir $BaseDir \n";
print data "unsetenv BetaModeString\n"; 
print data "unsetenv BetaColl\n";
print data "unsetenv lookMC\n";
print data "\n\n\n";

print "Creating the ChangeLogFilesName.script file \n";
open(change,">ChangeLogFilesName.script") or die;

#Now grep the files
my $file = "";
my $n = 0;

my @files = `grep -l xited *out`;
foreach $file (@files) {
    #At this point file is the name of the file + end-of-line
    my $f = `grep -lZ xited $file`;
    #now f is the file name without end-of-line

    #to change the name of the found files
    print change "mv $f $f-Error \n";
    #also doable this way? --> yes I think so...
    my @command = ("mv","$f","$f-Error");
    system(@command);

    #now decompose the name of the file to get the type-run-nb
    my $part = "";
    my $type = "";
    my $run = "";
    my $nb = "";
    
    $n=1;
    my @parts = split(/-/, $f);
    foreach $part (@parts) {
	if($n==1) { $type= "$part"; }
	elsif($n==2) { $run = "$part"; }
	elsif($n==3) {  
	    my $nn=0;
	    my @pps = split(/.out/, $part);
	    my $pp="";
	    foreach $pp (@pps) {
		if($nn==0) { $nb = "$pp"; }
		++$nn;
	    }
	}
	++$n;
    }

    #We now have the type-run-nb
    #Write it down to the appropriate file, so that it can be rerun
    if($type eq "BpBm" || $type eq "B0B0bar") { 
	print "Found $type $run $nb and wrote it to ../scripts/RerunJobs-genBB.script\n";
	print BB "echo Sending job $type-$run-$nb... \n";	
	print BB "setenv BetaCollFile $BaseDir/Skimmed_DATASETS/SP-$type-XSLBtoXulnuFilter-$run-R14-$nb.tcl \n";
	print BB "setenv RootFile  $BaseDir/ntuples/$type-$run-$nb.root \n";
	print BB "bsub -q long -o $BaseDir/logFiles/$type-$run-$nb.out ../bin/Linux24RH72_i386_gcc2953/BetaMiniApp ../XslReader/runReaderMC.tcl \n";
	print BB "\n";
    }
    elsif($type eq "uds" || $type eq "ccbar" || $type eq "tautau") { 
	print "Found $type $run $nb and wrote it to ../scripts/RerunJobs-Continuum.script\n";
	print cont "echo Sending job $type-$run-$nb... \n";	
	print cont "setenv BetaCollFile $BaseDir/Skimmed_DATASETS/SP-$type-XSLBtoXulnuFilter-$run-R14-$nb.tcl \n";
	print cont "setenv RootFile  $BaseDir/ntuples/$type-$run-$nb.root \n";
	print cont "bsub -q long -o $BaseDir/logFiles/$type-$run-$nb.out ../bin/Linux24RH72_i386_gcc2953/BetaMiniApp ../XslReader/runReaderMC.tcl \n";
	print cont "\n";
    }
    elsif($type eq "OffPeak" || $type eq "OnPeak") { 
	print "Found $type $run $nb and wrote it to ../scripts/RerunJobs-data.script\n";
	print data "echo Sending job $type-$run-$nb... \n";
	if(!($run eq "Run4")) { print data "setenv BetaCollFile $BaseDir/Skimmed_DATASETS/XSLBtoXulnuFilter-$run-$type-R14a-$nb.tcl \n"; }
	elsif(true){ print data "setenv BetaCollFile $BaseDir/Skimmed_DATASETS/XSLBtoXulnuFilter-$run-$type-R14-$nb.tcl \n"; }
	print data "setenv RootFile  $BaseDir/ntuples/$type-$run-$nb.root \n";
	print data "bsub -q long -o $BaseDir/logFiles/$type-$run-$nb.out ../bin/Linux24RH72_i386_gcc2953/BetaMiniApp ../XslReader/runReaderData.tcl \n";
	print data "\n";
    }
    else { print "Found $type $run $nb. This case is not taken care of by this script... Sorry!"; }

}

print BB "echo Voila \n\n";
close BB;

print cont "echo Voila \n\n";
close cont;

print data "echo Voila \n\n";
close data;


print "Voila! C'est termine! :-)) \n";
print "source ChangeLogFilesName.script \n";
print "will change the name of your problematic log files for you. You should do it before executing the other scripts... \n";




