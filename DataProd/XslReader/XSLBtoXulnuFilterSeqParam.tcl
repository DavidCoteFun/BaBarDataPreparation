
sourceFoundFile CompositionSequences/defineCompAlias.tcl

sourceFoundFile BetaMiniUser/btaMini.tcl
#sourceFoundFile PhysProdTools/LoadEventInfoSequence.tcl
sourceFoundFile BetaPid/PidMicroSequence.tcl
sourceFoundFile CompositionSequences/CompCharmlessProdSequence.tcl
sourceFoundFile CompositionSequences/CompPsiInitSequence.tcl
sourceFoundFile CompositionSequences/CompTrackRefineSequence.tcl
sourceFoundFile SimpleComposition/SmpCompositionSequence.tcl


#sequence disable CompMicroSequence
sequence disable CompKstarSequence

sequence enable PidMergedPi0MicroSequence
sequence enable PidMuonMicroSequence
sequence enable PidElectronMicroSequence
sequence enable PidKaonMicroSequence
sequence enable PidPionMicroSequence 
sequence enable CompEtaSequence
sequence enable CompTrackRefineSequence

module disable a1CToRho0Pi_Loose
module disable a1CToRho0Pi_Default
module disable a1CToRho0Pi_Tight
module disable  a1CToRho0PiHard_Tight

#module talk CompPi0ListMerger
#inputMergPi0 set AnotherName
#exit

sourceFoundFile ../XslReader/XSLBtoXulnuFilterSequence.tcl








