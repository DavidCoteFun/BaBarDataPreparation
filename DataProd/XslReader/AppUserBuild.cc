//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: AppUserBuild.cc,v 1.5 2005/11/09 17:48:03 denardo Exp $
//
// Description:
//	Class AppUserBuild implementation for the XslReader package.
//	This defines the link-time contents of the BetaApp application
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Bill Lockman
//      Akbar Mokhtarani
//
// Copyright Information:
//	Copyright (C) 2002		UC Santa Cruz, LBL
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include <iostream>
//-----------------------
// This Class's Header --
//-----------------------
#include "Framework/AppUserBuild.hh"
#include "Framework/AppFramework.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
// core sequences
#include "BetaMiniEnvSequences/BetaMiniEnvSequence.hh"
#include "BetaMiniSequences/BetaMiniSequence.hh"

// actions
#include "BetaMiniSequences/BetaMiniActions.hh"

// Minimilist set of physics sequences
#include "BetaMiniSequences/BetaMiniPhysicsSequence.hh"

// full PhysProdSequence
#include "PhysProdTools/PhysProdSequence.hh"

// QA module
#include "BetaMiniQA/BetaMiniQaSequence.hh"


//XSLBtoXulnu stuff:
#include "CompositionTools/CompMergeLists.hh"
#include "XslReader/XSLBtoXulnuFilterMini.hh"
#include "XslReader/XSLReader.hh"
#include "BetaPid/PidMergedPi0MicroSequence.hh"

//What's the difference between both?
#include "BetaPid/PidMicroSequence.hh"
#include "PidTools/PidSequence.hh"
#include "PhysProdTools/LoadEventInfoSequence.hh"

//Seq
#include "CompositionSequences/CompCharmlessProdSequence.hh"
#include "CompositionSequences/CompPsiInitSequence.hh"
#include "CompositionSequences/CompTrackRefineSequence.hh"


//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//----------------
// Constructors --
//----------------

AppUserBuild::AppUserBuild( AppFramework* theFramework )
    : AppBuild( theFramework )
{
  // core sequence
  BetaMiniEnvSequence(this);
  BetaMiniSequence(this);

  // actions
  BetaMiniActions(theFramework);

  // QA sequence
  BetaMiniQaSequence(this);

  // minimilist physics sequence
  BetaMiniPhysicsSequence(this);


  //XSLBtoXulnu Filter and utilities
  add(new XSLBtoXulnuFilterMini("XSLBtoXulnuFilterMini", "filter module for exlusive B->Xulnu decays"));
  add(new XSLReader("XSLReader", "analysis module reading the output of XSLBtoXulnuFilter"));

  //PID with minimal physics sequence
  PidSequence( this );
  PidMicroSequence( this );
  LoadEventInfoSequence( this );

  //Pi0's
  PidMergedPi0MicroSequence(this); 

  //dispatching the new lists in CompositionSequences
  CompCharmlessProdSequence(this);  //rho's and omega
  CompPsiInitSequence(this);  //BremRecovery
  CompTrackRefineSequence(this);  //Track refine sequence

}

//--------------
// Destructor --
//--------------

AppUserBuild::~AppUserBuild( )
{
}

