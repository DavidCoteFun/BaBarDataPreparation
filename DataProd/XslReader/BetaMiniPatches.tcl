# $Id: BetaMiniPatches.tcl,v 1.1 2004/10/05 23:41:39 cote Exp $
#------------------------------------------------------------------------------

#if [ info exists ConfigPatch ] {
#    if [ string match $ConfigPatch "MC" ] {
## patch needed to run on SP4 MC
#	talkto DchBuildEnv {
#	    dedxCalibBkgStrategy set f
#	    puts "BetaMiniPatches: Setting DchBuildEnv::dedxCalibBkgStrategy to false"
#	}
#    }
#}

## a patch needed to run on 10 series schema
#talkto SvtBuildEnv {
#    SvtWaferDistortionsFile set SvtGeom/NullDistortion.dat
#}

ErrMsg trace " completed OK"
