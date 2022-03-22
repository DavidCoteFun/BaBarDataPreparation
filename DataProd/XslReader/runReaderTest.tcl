#------------------------------------------------------------------------------
# $Id: runReaderTest.tcl,v 1.1 2005/07/16 05:17:34 cote Exp $
# Sample MyMiniAnalysis.tcl file
#------------------------------------------------------------------------------
# always source the error logger early in your main tcl script
sourceFoundFile ErrLogger/ErrLog.tcl
sourceFoundFile FrameScripts/talkto.tcl

sourceFoundFile CompositionSequences/defineCompAlias.tcl

# set the error logging level to 'warning'.  If you encounter a configuration error
# you can get more information using 'trace'
ErrLoggingLevel warning


## allowed values of BetaMiniReadPersistence are (currently) "Kan", "Bdb"
##
set BetaMiniReadPersistence Kan
#set BetaMiniReadPersistence Bdb

## allowed (non-expert) values of levelOfDetail are "micro", "cache", "extend" or "refit"
##
set levelOfDetail "cache"

## allowed values of ConfigPatch are "Run1", "Run2" or "MC".  This MUST be set consistent
## with your input data type or you will get INCONSISTENT OR INCORRECT RESULTS
##

FwkCfgVar ConfigMC $env(DumpMCInfo)

if { $ConfigMC == "true"}  {
    set ConfigPatch "MC"
} else {
    set ConfigPatch "Run2"
}


##
##  You can enter input collections two ways: either append them to a list, or
##  explicitly enter them in the input module.  Do one or the other, BUT NOT BOTH.
##  If inputList is set before executing btaMini.tcl, that will automatically
##  add the collections to the appropriate input module, otherwise make sure you
##  talk to the right one.
##
## lappend inputList collection1 collection2 ...

if [ info exists env(BetaCollFile) ] {
    source $env(BetaCollFile)
} else {
    echo BetaCollFile not set, not configuring an input file
}


if [ info exists env(BetaColl) ] {
    lappend inputList $env(BetaColl)
} else {
    echo BetaColl not set, not configuring an inputList
}


## choose the flavor of ntuple to write (hbook or root) and the file name
##
set BetaMiniTuple "root"
set histFileName $env(RootFile)

## create Everything path and add core sequences to it.  btaMiniPhysics is the same as btaMini, just
## appending a few standard list generating modules.  for reading data with stored composites, you may
## have a conflict running btaMiniPhyscs.tcl
##
## You can also run (most of) the PhysProdSequence, complete with its 3 gamma conversion
## finders, etc.  Consider disabling the portion of this sequence that
## you do not need to save yourself some time.  The BetaLumiSequence
## and TagProd sequences are left off, as they otherwise cause problems.
##


#Setting up the path
sourceFoundFile BetaMiniUser/btaMini.tcl
#sourceFoundFile BetaMiniUser/btaMiniPhysics.tcl
#sourceFoundFile BetaMiniUser/btaMiniPhysProdSequence.tcl
#path append Everything XSLReader


sourceFoundFile ../XslReader/XSLReaderSequence.tcl
path append Everything XSLReaderSequence


#Reading back UsrData
sourceFoundFile UsrTools/UsrDataProcs.tcl
enableReadUsrData
readEventUsrData XSLBtoXulnuEventData


if { $ConfigMC == "true"}  {
    
    ##Setting up PID tables and killing/tweaking/weighting
#    pidCfg_dataset cdb 
    pidCfg_dataset ascii
    
    #the ones set to weight are those used by XSLBtoXulnuFilter
    pidCfg_mode tweak VeryLooseElectronMicroSelection
    pidCfg_mode tweak LooseElectronMicroSelection
    pidCfg_mode tweak TightElectronMicroSelection
    pidCfg_mode tweak VeryTightElectronMicroSelection
    pidCfg_mode weight PidLHElectronSelector
    
    pidCfg_mode tweak TightMuonMicroSelection
    pidCfg_mode tweak VeryTightMuonMicroSelection
    pidCfg_mode tweak NNVeryLooseMuonSelectionFakeRate
    pidCfg_mode tweak NNVeryLooseMuonSelection
    pidCfg_mode tweak NNLooseMuonSelectionFakeRate
    pidCfg_mode tweak NNLooseMuonSelection
    pidCfg_mode tweak NNTightMuonSelectionFakeRate
    pidCfg_mode weight NNTightMuonSelection
    pidCfg_mode tweak NNVeryTightMuonSelectionFakeRate
    pidCfg_mode weight NNVeryTightMuonSelection
    
    pidCfg_mode weight LooseLHPionMicroSelection
    pidCfg_mode tweak TightLHPionMicroSelection
    pidCfg_mode tweak VeryTightLHPionMicroSelection
    
    pidCfg_mode tweak LooseLHKaonMicroSelection
    pidCfg_mode tweak TightLHKaonMicroSelection
    pidCfg_mode tweak VeryTightLHKaonMicroSelection

    pidCfg_mode tweak LooseLHProtonSelection
    pidCfg_mode tweak TightLHProtonSelection
    pidCfg_mode tweak VeryTightLHProtonSelection

    ##Neutrals corrections
    talkto EmcNeutCorrLoader {
	correctionOn set true
	endcapShift  set true }
    
}


##
##  If your job has a tag-level filter, here is how you should run it
##  so as to avoid wasting time reading the mini when the tag filter fails
##  Here's a simple example that restricts to just multi-hadron events
##  on Kan input

module clone TagFilterByName TagBGFMultiHadron
module talk TagBGFMultiHadron
  andList set BGFMultiHadron
  assertIfMissing set true
exit
sequence append BetaMiniReadSequence -a KanEventUpdateTag TagBGFMultiHadron

mod talk EvtCounter
verbose set true
#skipStart set $env(startSkip)
#skipStop  set $env(endSkip)
#badEvents set $env(BadEvts)
printFreq set 1
print set true
exit

mod talk XSLReader
MyVerbose set false
lookMCTruth set $env(DumpMCInfo)
readingMode set Run4
signalModeString set $env(BetaModeString)
AnalyzePilnu set $env(AnalPilnu)
AnalyzePi0lnu set $env(AnalPi0lnu)
AnalyzeEtalnu set $env(AnalEtalnu)
AnalyzeEtaplnu set $env(AnalEtaplnu)
AnalyzeRhoClnu set $env(AnalRhoClnu)
AnalyzeRho0lnu set $env(AnalRho0lnu)
AnalyzeOmegalnu set $env(AnalOmegalnu)
AnalyzeGammalnu set false
exit


path list
ev begin -nev 1000

ErrMsg trace "completed OK"
exit







