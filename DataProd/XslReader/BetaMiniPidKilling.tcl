# $Id: BetaMiniPidKilling.tcl,v 1.1 2004/10/05 23:41:39 cote Exp $
# 
# TCL script to run the PID-Killing procedure with the "Mini".
# 
#------------------------------------------------------------------------------


# Check if Killing enviorement variables exist
if {![info exists BtaPidKillMuon] &&
    ![info exists BtaPidKillElectron] &&
    ![info exists BtaPidKillKaon] &&
    ![info exists BtaPidKillPion] &&
    ![info exists BtaPidKillProton] &&
    ![info exists BtaPidKilling] } {

ErrMsg fatal "No PID-Killing configuration variables set in this job."

}

# Run Killing
global PidKill
set PidKill no
if { [info exists BtaPidKilling] } {
    set PidKill $BtaPidKilling
}
if { $PidKill=="yes" } {
    ErrMsg warning "PID-Killing Enable for All Particles."
}
global PidKillSingle
set PidKillSingle no

sequence create PidKillingMicroSequence
#sequence insert EventSource -b PhysBetaMicroSequence BetaLoadPidKillingMicroSequence
sequence append PidProdSequence PidKillingMicroSequence
sequence append PidKillingMicroSequence BtaMicroPidKilling

sourceFoundFile BetaMicro/particleKill.tcl

#########
# Muons #
#########

if [ info exists BtaPidKillMuon ] {
    if { $BtaPidKillMuon=="yes" } {
	set PidKillSingle yes
	ErrMsg warning "PID-Killing Enable for Muons."
    } else {
	set PidKillSingle no
    }
} else {
    set PidKillSingle no
}
foreach type { MinimumIonizing VeryLoose Loose Tight VeryTight } {
    lappend modules ${type}MuonMicroSelection
    lappend lists muMicro${type}
}

if { $PidKill=="yes" || $PidKillSingle=="yes" } {
    particleKill Muon mu $modules $lists
    unset modules
    unset lists
}

#############
# Electrons #
#############
if [ info exists BtaPidKillElectron ] {
    if { $BtaPidKillElectron=="yes" } {
	set PidKillSingle yes
	ErrMsg warning "PID-Killing Enable for Electrons."
    } else {
	set PidKillSingle no
    }
} else {
    set PidKillSingle no
}
foreach type { NoCal VeryLoose Loose Tight VeryTight } {
    lappend modules ${type}ElectronMicroSelection
    lappend lists eMicro${type}
}
# LH selector naming convention is a little different
lappend modules PidLHElectronSelector
lappend lists PidLHElectrons

if { $PidKill=="yes" || $PidKillSingle=="yes" } {
    particleKill Electron e $modules $lists
    unset modules
    unset lists
}

#########
# Kaons #
#########
if [ info exists BtaPidKillKaon ] {
    if { $BtaPidKillKaon=="yes" } {
	set PidKillSingle yes
	ErrMsg warning "PID-Killing Enable for Kaons."
    } else {
	set PidKillSingle no
    }
} else {
    set PidKillSingle no
}
foreach type { VeryLoose Loose Tight VeryTight NotPion } {
    # micro selectors
    lappend modules ${type}KaonMicroSelection
    lappend lists KMicro${type}
    # NN selectors
    lappend modules ${type}NNKaonMicroSelection
    lappend lists KNN${type}
    # LH selectors
    lappend modules ${type}LHKaonMicroSelection
    lappend lists KLH${type}
}

if { $PidKill=="yes" || $PidKillSingle=="yes" } {
    particleKill Kaon K $modules $lists
    unset modules 
    unset lists 
}

###########
# Protons #
###########
if [ info exists BtaPidKillProton ] {
    if { $BtaPidKillProton=="yes" } {
	set PidKillSingle yes
	ErrMsg warning "PID-Killing Enable for Protons."
    } else {
	set PidKillSingle no
    }
} else {
    set PidKillSingle no
}
# LH selection
foreach type { Loose VeryLoose Tight VeryTight } {
    lappend modules ${type}LhProtonSelection
    lappend lists pLH${type}
}

if { $PidKill=="yes" || $PidKillSingle=="yes" } {
    particleKill Proton p $modules $lists
    unset modules 
    unset lists 
}
# GRL selectors (aka Micro selectors) deprecated
#  because EMC/DRC consistencies aren't stored correctly
#  # Grl selectors (make "pMicro" lists)
#  foreach type { Loose Default Tight } {
#      lappend modules Grl${type}ProtonSelection
#      lappend lists pMicro${type}
#  }

#########
# Pions #
#########
if [ info exists BtaPidKillPion ] {
    if { $BtaPidKillPion=="yes" } {
	set PidKillSingle yes
	ErrMsg warning "PID-Killing Enable for Pions."
    } else {
	set PidKillSingle no
    }
} else {
    set PidKillSingle no
}
# Roy selectors
foreach type { Loose NotKaon } {
    lappend modules PidRoyPionSelection${type}
    lappend lists piRoy${type}
}
# Lh selection
foreach type { VeryLoose Loose Tight VeryTight } {
    lappend modules ${type}LHPionMicroSelection
    lappend lists piLH${type}
}

if { $PidKill=="yes" || $PidKillSingle=="yes" } {
    particleKill Pion pi $modules $lists
    unset modules 
    unset lists 
}
