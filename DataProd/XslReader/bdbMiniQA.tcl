# $Id: bdbMiniQA.tcl,v 1.1 2004/10/05 23:41:44 cote Exp $
source PARENT/FrameScripts/sourceFoundFile.tcl
sourceFoundFile FrameScripts/setProduction.tcl
sourceFoundFile FrameScripts/talkto.tcl

puts "##### bdbMiniQA.tcl is deprecated. Please source BetaMiniQA/btaMiniQA.tcl instead"
sourceFoundFile BetaMiniQA/btaMiniQA.tcl
puts "bdbMiniQA.tcl completed with WARNINGS"