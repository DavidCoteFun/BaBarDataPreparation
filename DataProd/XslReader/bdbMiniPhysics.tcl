# $Id: bdbMiniPhysics.tcl,v 1.1 2004/10/05 23:41:44 cote Exp $
source PARENT/FrameScripts/sourceFoundFile.tcl
sourceFoundFile FrameScripts/setProduction.tcl
sourceFoundFile FrameScripts/talkto.tcl

puts "##### bdbMiniPhysics.tcl is deprecated. Please source btaMiniPhysics.tcl instead"
sourceFoundFile BetaMiniUser/btaMiniPhysics.tcl
puts "BetaMiniUser/bdbMiniPhysics.tcl completed with WARNINGS"
