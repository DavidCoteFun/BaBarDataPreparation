# $Id: btaMiniPhysProdSequence.tcl,v 1.1 2004/10/05 23:41:44 cote Exp $
# btaMiniPhysics.tcl:
# added BetaMiniPhysicsSequence on top of btaMini.tcl
#------------------------------------------------------------------------------

sourceFoundFile BetaMiniUser/btaMini.tcl

## add the kitchen sink
##

global fullMiniPhysSequence
set fullMiniPhysSequence 1

sourceFoundFile BetaMiniSequences/BetaMiniPhysicsSequence.tcl
path append Everything BetaMiniPhysicsSequence

ErrMsg trace  " completed OK"
