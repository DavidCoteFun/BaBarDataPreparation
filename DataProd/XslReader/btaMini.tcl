# $Id: btaMini.tcl,v 1.1 2004/10/05 23:41:44 cote Exp $
# btaMini.tcl:
# TCL script to run BetaMiniApp executable to read the "Mini" contents
# from the event store and perform user analysis
#------------------------------------------------------------------------------
## source the FrameScripts stuff for others to use
##
sourceFoundFile FrameScripts/setProduction.tcl
sourceFoundFile FrameScripts/talkto.tcl
sourceFoundFile ErrLogger/ErrLog.tcl

# put together all the pieces in a path
path delete AllPath
path create Everything 

## insert core sequence
##
sourceFoundFile BetaMiniSequences/BetaMiniSequence.tcl 
path append Everything BetaMiniSequence

# look for a patches file in the current directory
set patches BetaMiniPatches.tcl
if [ file exists ./$patches ] {
    set patchesfile ./$patches
} 
if [ info exists patchesfile ] {
    ErrMsg warning  "sourcing $patchesfile"
    source $patchesfile
}

ErrMsg routine  "completed OK"
