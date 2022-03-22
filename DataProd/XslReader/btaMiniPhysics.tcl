# $Id: btaMiniPhysics.tcl,v 1.1 2004/10/05 23:41:44 cote Exp $
# btaMiniPhysics.tcl:
# added BetaMiniPhysicsSequence on top of btaMini.tcl
#------------------------------------------------------------------------------

sourceFoundFile BetaMiniUser/btaMini.tcl

## add the minimal physics list
##
sourceFoundFile BetaMiniSequences/BetaMiniPhysicsSequence.tcl
path append Everything BetaMiniPhysicsSequence

ErrMsg trace "completed OK"
