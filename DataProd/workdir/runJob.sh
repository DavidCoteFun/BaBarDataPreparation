echo Producing ntuples for failed jobs

unsetenv BetaColl
setenv BetaModeString data 
setenv DumpMCInfo false
setenv AnalPilnu true 
setenv AnalPi0lnu false 
setenv AnalEtalnu false 
setenv AnalEtaplnu false 
setenv AnalRhoClnu false 
setenv AnalRho0lnu false 
setenv AnalOmegalnu false 

echo Sending job OnPeak-Run5-R18b-v03-30... 
setenv BetaCollFile /nfs/slac/g/udem/u3/users/ProdV10/Skimmed_DATASETS_v03/XSLBtoXulnuFilter-Run5-OnPeak-R18b-v03-30.tcl 
setenv RootFile  /nfs/slac/g/udem/u3/users/ProdV10/pilnu/ntuples_v03/OnPeak-Run5-R18b-v03-30.root 
bsub -q kanga -o /nfs/slac/g/udem/u3/users/ProdV10/pilnu/logFiles_v03/OnPeak-Run5-R18b-v03-30.out ../bin/Linux24SL3_i386_gcc323-noOptimize-Debug/BetaMiniApp ../XslReader/runReader.tcl 

echo voila


