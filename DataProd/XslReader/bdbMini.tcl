# $Id: bdbMini.tcl,v 1.1 2004/10/05 23:41:44 cote Exp $
source PARENT/FrameScripts/sourceFoundFile.tcl
sourceFoundFile FrameScripts/setProduction.tcl
sourceFoundFile FrameScripts/talkto.tcl

puts "##### bdbMini.tcl is deprecated. Please source btaMini.tcl instead"
sourceFoundFile BetaMiniUser/btaMini.tcl
puts "BetaMiniUser/bdbMini.tcl completed with WARNINGS"
