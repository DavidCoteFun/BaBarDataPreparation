#------------------------------------------------------------------------------
# $Id: runFilterDATA.tcl,v 1.1 2004/10/05 23:41:44 cote Exp $
# Sample MyMiniAnalysis.tcl file
#------------------------------------------------------------------------------
# always source the error logger early in your main tcl script
sourceFoundFile ErrLogger/ErrLog.tcl
source PARENT/FrameScripts/sourceFoundFile.tcl
sourceFoundFile FrameScripts/talkto.tcl


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
set ConfigPatch "Run2"

##
##  You can enter input collections two ways: either append them to a list, or
##  explicitly enter them in the input module.  Do one or the other, BUT NOT BOTH.
##  If inputList is set before executing btaMini.tcl, that will automatically
##  add the collections to the appropriate input module, otherwise make sure you
##  talk to the right one.
##
## lappend inputList collection1 collection2 ...

if [ info exists env(CollFile) ] {
    source $env(CollFile)
} else {
    echo CollFile not set, not configuring an input file
}


if [ info exists env(Coll) ] {
    lappend inputList $env(Coll)
} else {
    echo Coll not set, not configuring an inputList
}


##
##  OR THE FOLLOWING (choose the correct one based on persistence)
##
## talkto BdbEventInput {
## talkto KanEventInput {
##    input add collection1
##    input add collection2
##    ...
## }

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
#sourceFoundFile BetaMiniUser/btaMini.tcl
#sourceFoundFile BetaMiniUser/btaMiniPhysics.tcl

#Uncomment this to run the "kitchen sink"...
#sourceFoundFile BetaMiniUser/btaMiniPhysProdSequence.tcl
#path append Everything XSLBtoXulnuFilterMini
#path append Everything XSLBtoXulnuReader
 
#Uncomment this to run XSLBtoXulnuFilterSequence...
sourceFoundFile ../XslReader/XSLBtoXulnuFilterSeqParam.tcl
path append Everything XSLBtoXulnuFilterSequence

##
##  If your job has a tag-level filter, here is how you should run it
##  so as to avoid wasting time reading the mini when the tag filter fails
##  Here's a simple example that restricts to just multi-hadron events
##  on Kan input

#module clone TagFilterByName TagBGFMultiHadron
#module talk TagBGFMultiHadron
#  andList set BGFMultiHadron
#  assertIfMissing set true
#exit
#sequence append BetaMiniReadSequence -a KanEventUpdateTag TagBGFMultiHadron

mod talk EvtCounter
verbose set true
#skipStart set 10001
#skipStop  set 20000
printFreq set 1
print set true
exit


## configure MyMiniAnalysis module
##
mod talk XSLBtoXulnuFilterMini
Verbose set true
exit

mod talk XSLBtoXulnuReader
MyVerbose set false
AnalyzeGamma set true
readFromFilterMini set true
exit


path list
ev begin -nev 50
#ev begin -nev $env(NbEv)

ErrMsg trace "completed OK"
exit










