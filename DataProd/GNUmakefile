# top level GNUmakefile for a release
# $Id: GNUmakefile.main,v 1.551 2005/12/05 09:51:57 terryh Exp $
#
# Note: See standard.mk for more variable definition/descriptions.
#
#___________________________________________________________________________
#++ Major Targets:
#
#    all	equivalent to 'lib', followed by 'bin'
#    bin        build binaries and scripts defined in BINS and BINSCRIPTS.  
#    bootstrap  build basic packages that need to be built in special order.
#    clean	delete files in working directories for all architectures.
#    cold      used when making a production release; bootstrap + all.
#    installdirs   create needed directories.
#    lib        build library.
#    libext     build library that can be loaded independantly of build.
#
#++ Minor targets:
#
#    allbin     build executables defined in BINS + EXTRABINS + BINSCRIPTS +EXTRABINSCRIPTS.
#    binscripts build scripts defined in BINSCRIPTS.
#    cleanarch	delete files in working directories for current $BFARCH.
#    cleanbin   clean binary only.
#    cleanlib   clean library only.
#    extrabin   build extra binaries defined in EXTRABINS.
#    extrabinscripts   build extra scripts defined in EXTRABINSCRIPTS.
#    pkg.<foo>  build individual binary with 'gmake pkg.foo'.
#    <foo>      'gmake <foo>' in all subdirectories. <foo> not containing a .)
#    siteinstall install scripts for remote site. Normally run binscript.
#
# Typically, a newly created release within the BaBar system
#  needs a gmake installdirs, followed by one or more gmakes.
#
# One can also do particular directories, like 'gmake stdhep.lib'
#  or 'gmake HepTuple.bin',  or a complete phase like 'gmake lib'
#
# The include path now starts with tmp/$BFARCH, then <>/$BFCURRENT/tmp
# so that locally created include files will be found 
# preferentially
#
#___________________________________________________________________________
#++ Options that can be set in BFARCH or ROPT. e.g. "gmake ROPT=-Time-Optimize"
#   Use 'no' to turn off an option; e.g. ROPT=-noDebug.
#
#++ Major ROPT options:
#    -Debug	turn on debug option. e.g. add '-g' to CXXFLAGS.
#    -NOBIN	do not build binary except in bootstrap
#    -Optimize  turn on optimization. e.g. add '-O' to CXXFLAGS.
#    -Production  production release. Same as PRODUCTION_RELEASE=1.
#    -Prunebin    prune binary based on keep_xxx list.
#    -Shared    build shared library. =-Shared-Lshared-Lstatic.
#    -Verbose   Verbose mode, same as VERBOSE=1.
#
#++ Minor ROPT options:
#    -Fastbuild collect .o to libfoo.a at end of lib target to save I/O.
#    -Ldlink2   uses ldlink to find shared SYSLIB, turn off "-rpath/-R".
#    -Linkmap   generates link map in bin/$BFARCH/<binary>.map.
#    -Lshared   put shared library in link library search path.
#    -Lstatic   put static library in link library search path, default.
#    -native	use native compiler.
#    -NOInstalldirs  skip installdirs
#    -Objy	use Objectivity.
#    -Profile	profile build. Add -pg to CXXFLAGS etc.
#    -Showdir   show detail directory level information pkg.dir.target"
#    -Showdep   show dependents of binaries (based on -l<name> in LOADLIBES)
#    -Skipdfile skip 'include .d files'
#    -Skiplib   skip 'lib' step. Use production library.
#    -SkipSlaclog skip 'slaclog' when it hangs.
#    -Static    build static library, default. Exclusive with -Shared.
#    -Time	display CPU timing information.
#
#___________________________________________________________________________
#++ Options used in 'gmake clean' and 'gmake cleanarch'
#
#    CLEANLIST  full path of directories to be cleaned.
#		e.g CLEANLIST='lib tmp/SunOS5'  
#    NOCLEAN    top directories not to be cleaned. e.g. NOCLEAN='lib results'
#
# Note:
#    CLEANLIST must be the full path while NOCLEAN is the top-level directory. 
#
#___________________________________________________________________________
#++ Variables used to define general directores.
#
#    TOPDIR	top level directory of this production/test release.
#    srcdir	current directory in packages/sub-directories.
#    PACKAGE    name of package
#
#___________________________________________________________________________
#++ Centrally controlled list (Major)
#
# These lists are now contained in the babar_specials.mk file in SRT
#
# *) build list:
#    binonly_sunos5	build 'bin' for this package on SunOS5 only.
#    nodebug_pkg	packages to be built with -noDebug
#    debug_pkg		packages to be built with -Debug
#    noautomatic_pkg	packages to be built with -fno-automatic
#    only_sunos5 build 'all' for this package on SunOS5 only. For online build.
#    skip_optimize	do not optimize these packages in -Optimize build.
#
# *) keep/nokeep list:
#    keep_bin		binaries to be kept in a -Prunebin build	
#    keep_bin_opt	binaries to be kept in a -Prunebin-Optimize build
#    keep_pkg		packages to be kept in a -Prunebin build	
#    keep_pkg_opt	packages to be kept in a -Prunebin-Optimize build
#    nokeep_bin		binaries that should be pruned in all cases
#    nokeep_bin_opt	binaries that should be pruned in an -Optimize build
#    nokeep_bin_noopt	binaries that should be pruned in an -noOptimize build
#
#   Centrally controlled list (Minor)
#    protect_dbdir	database directory to be protected	
#    skip_bin		skip bin stage
#    skip_fastbuild	ignore -Fastbuild option for these packages
#    skip_include	skip include stage
#    skip_include_d	skip include *.d files. e.g. bootstrap, clean 
#    skip_lib		skip lib stage
#    skip_nonstd	do not use -fnonstd flag to dump upon FPE.
#    skip_test		skip test stage
# 
#___________________________________________________________________________
#   

#___________________________________________________________________________
# Debugging Options:
#
# VERBOSE=yes or 1  Dump of the rules and other verbose output.
# DEBUG=1 : debug ROPT, OPTIONS.
# DEBUG02=1 : debug schema information.
# DEBUG03=1 : debug subdirectory level, turn on 'set -x'. 
# DEBUG04=1 : debug sub-make, MAKEOVERRIDES. See standard.mk.
# DEBUG05=1 : debug shared library. See standard.mk.
# DEBUG06=1 : debug ALLLIB of shared library, careful! long! See standard.mk.
# DEBUG07=1 : debug ROPT=-Prunebin option.
# DEBUG08=1 : debug SRCDIR, etc
# DEBUG10=1 : debug include .d files. See standard.mk.
# DEBUG11=1 : debug speed by monitoring PACKAGES, etc.
# DEBUG12=1 : debug checking links during installdirs
# DEBUG13=1 : debug display BIN1, LIBCCFILE1, etc. See standard.mk.
# DEBUG14=1 : debug corba.
# DEBUG15=1 : debug 'clean'.
# DEBUG16=1 : debug slaclog. Won't call slaclog. Use 'logopt=-debug' to debug.
# DEBUG17=1 : debug binonly-sunos5.
# DEBUG18=1 : debug nodebug_pkg and debug_pkg
# DEBUG19=1 : test nodebug_pkg, define this list for test release
# DEBUG20=1 : test clean target's NOCLEAN and CLEANLIST option
# DEBUG21=1 : debug BIN targets, display BINNAME and MAKEGOAL.
# DEBUG23=1 : debug message of installdirs
# DEBUG24=1 : debug compiler checking.
# DEBUG25=1 : debug checking of main().
# DEBUG26=1 : debug ldlink, major messages only.
# DEBUG27=1 : debug ldlink, detail messages.
# DEBUG28=1 : run ldlink by running "<pkg>.ldlink" for each package. Avoid it!
# DEBUG29=1 : trace entrance to GNUmakefile, arch_spec.mk, and standard.mk.
# DEBUG30=1 : turn on informational messages for ldlink, normally suppressed.
# DEBUG31=1 : debug message for srtbinlist.
# DEBUG32=1 : debug message for rootlib (multi-package ROOT library) 
# DEBUG33=1 : debug variables for arch_spec_flex.mk.
# DEBUG34=1 : debug mkdir tmp and shtmp

# BobJ Oct 94 - debugging version
# BobJ Mar 95 - reworked to add general target passdown to subdirs
# BobJ June 95 - removed FrameExample default link step (see TODO)
# BobJ Dec 95 - changed logic for finding $PACKAGES
#          to take (any and only) subdirectories containing a GNUmakefile
#          instead of explicit exception list.  Removed the explicit
#          FrameExample.bin defeat, which should have been limited to
#          LLNL only (other HP sites have bothered to install Tcl)
# T. Wenaus Nov 95 - rule printing and 'entering dir' printing turned off
#	 by default. Can be activated via VERBOSE=yes
# B. Lockman Mar 96 - don't do dbx if gnbase is r/o
# M.Ronan  July 96 - Add commissioning detectors to include and bin lists.
# B. Lockman Sep 96 - replaced dbx target in all with debugger which now
#	builds .dbxinit and .gdbinit directory files for dbx and gdb 
#	debuggers, respectively.
# H. Park  Jan 97 - Add printing lines in log files  for automations
# M. Kelsey  May 98 -- Replace use of shell to get file lists with use of
#	GMake's own $(wildcard).
# M.Kunze July 98 - Add support for shared libs
# S. Metzler July 98 - add DhpTuple DistObjTree to bootstrap
# S. Metzler April 99 - add HepScaler to bootstrap
# gpdf June 99 - add decodeOdfTime to keep_bin
# A. Samuel September 99 - add OepFastMonRemoteApp to keep_bins
# G. Zioulas September 99 - add OdcNtupleMaker to keep_bin
# E. Chen October 99 - add DhpMonControl to SunOS5 only
# S. Gowdy Novemver 99 - Add oocheckdbexist to keep_bin
# J. Becla November 99 - Add oodeletecont BdbBrowseFullDbs BdbConvertDbIdToName
#	to keep_bin
#

#==========================================================================
#++ SRT Initialization for each invocation of gmake, sub-gmake 

ifdef DEBUG29
   ERROR := $(shell echo "-> GNUmakefile $(MAKECMDGOALS)" >&2)
endif

SHELL = /bin/sh

ifeq (,$(TOPDIR))
ifeq ($(BFSITE),slac)
    export TOPDIR := $(shell a2nfs -noe -t `/bin/pwd` 2>/dev/null)
endif
ifeq (,$(TOPDIR))
    export TOPDIR := $(shell /bin/pwd)
endif
endif

export PWD := $(TOPDIR)
export MAKEBASE := $(TOPDIR)

# which SoftRelTools add-on is being used?
ifneq (,$(wildcard $(TOPDIR)/SoftRelTools))
  SRTPATH = $(TOPDIR)/SoftRelTools
else
ifneq (,$(wildcard $(BFOVERRIDE)/$(BFCURRENT)/SoftRelTools))
  SRTPATH = $(BFOVERRIDE)/$(BFCURRENT)/SoftRelTools
else
  SRTPATH = $(BFDIST)/releases/$(BFCURRENT)/SoftRelTools
endif
endif
export SRTPATH

# which full SoftRelTools is being used?
ifneq (,$(wildcard $(TOPDIR)/SoftRelTools))
  SRTPATH2 = $(TOPDIR)/SoftRelTools
else
  SRTPATH2 = $(BFDIST)/releases/$(BFCURRENT)/SoftRelTools
endif
export SRTPATH2


ifneq (,$(findstring Verbose,$(ROPT)))
  VERBOSE=1
endif
ifneq ($(VERBOSE),)
  MAKEFLAGS =
else
  MAKEFLAGS = --no-print-directory -s
endif

ifeq (,$(findstring SunOS,$(BFARCH)))
  ETIME=/usr/bin/time
else
  ETIME=/bin/time
endif

#__________________________________________________________________________
#++ findmk: subroutine to search for a mk/mk+ file
#   please provide 2 gmake macros: PKG: package MK: files to be included

define findmk
	if  [ -f $(BFOVERRIDE)/$(BFCURRENT)/$(PKG)/$(MK) ]; then \
	   echo $(BFOVERRIDE)/$(BFCURRENT)/$(PKG)/$(MK) ; \
	elif [ -f $(TOPDIR)/SiteConfig/$(BFSITE)/$(PKG)/$(MK) ]; then \
	   echo $(TOPDIR)/SiteConfig/$(BFSITE)/$(PKG)/$(MK) ; \
	elif [ -f $(BFDIST)/releases/$(BFCURRENT)/SiteConfig/$(BFSITE)/$(PKG)/$(MK) ]; then \
	   echo $(BFDIST)/releases/$(BFCURRENT)/SiteConfig/$(BFSITE)/$(PKG)/$(MK); \
	elif [ -f $(TOPDIR)/$(PKG)/$(MK) ]; then \
	   echo $(TOPDIR)/$(PKG)/$(MK); \
	elif [ -f $(BFDIST)/releases/$(BFCURRENT)/include/$(PKG)/$(MK) ]; then \
	   echo $(BFDIST)/releases/$(BFCURRENT)/include/$(PKG)/$(MK); \
	else \
	   echo "$(MK)"; \
	fi;
endef

#__________________________________________________________________________
#++ function to log monitoring information via slaclog
#   syntax: "slaclog bfbuild <start/end> <rel> <production/test> <goal> <arg>"
#   slaclog1: start record
#   slaclog2: end record

ifneq (,$(findstring noSkipSlaclog,$(ROPT)))
ifeq (,$(NOTFIRST))
override debug=
logopt=
ifdef DEBUG16
    override debug=echo
endif

#++ [slaclog1]: start record
define slaclog1
	if [ -r /usr/local/bin/slaclog ]; then \
	    $(debug) /usr/local/bin/slaclog $(logopt) bfbuild "start `cat .current` $(buildtype) $@ $(TOPDIR) $(MAKEOVERRIDES)"; \
	fi
endef

#++ [slaclog2]: end record
define slaclog2
	if [ -r /usr/local/bin/slaclog ]; then \
	    $(debug) /usr/local/bin/slaclog $(logopt) bfbuild "end `cat .current` $(buildtype) $@ $(TOPDIR)"; \
	fi
endef
endif   ## ~NOTFIRST
endif   ## ~SkipSlaclog

#==========================================================================
#++ Process build options. One time only section. This section is processed
#   only once for test release, twice for production/nightly release.

ifeq (,$(NOTFIRST))
    PKG = SoftRelTools
    MK=OS_check.mk
    include $(shell $(findmk))

     export BFARCHBASE := $(firstword $(subst -, ,$(BFARCH)))

#__________________________________________________________________________
#++ setting DEFOPT, default for all users

    override DEFOPT := -Debug-native-Objy-Optimize-Fastbuild-Ldlink2-SkipSlaclog

# Check ldlink2 in production  and give warning   (will be obsolete)
#ifeq (,$(PRODUCTION_RELEASE)$(PRODUCTION_TEST))
#  ifeq (,$(findstring ldlink,$(MAKECMDGOALS)))
#  ifeq (1,$(shell if [ "`ls ./shlib/$(BFARCH)/ldlink/*_pkg* $(BFDIST)/releases/$(BFCURRENT)/shlib/$(BFARCH)/ldlink/*_pkgid_* $(BFDIST)/releases/$(BFCURRENT)/BASE/shlib/$(BFARCH)/ldlink/*_pkgid_* 2>/dev/null`" = "" ]; then echo 1; fi))
#      error := $(shell echo "** Warning: no ldlink2 in production, please run 'gmake ldlink'." >&2)
#  endif   # no ldlink2
#  endif   # not ldlink
#endif     # not production

ifneq (,$(PRODUCTION_RELEASE)$(PRODUCTION_TEST))
    override buildtype := production
    override DEBUG30 := 1
else
    override buildtype := test
endif

#__________________________________________________________________________
#++ platform specific section, for all users

ifneq (,$(findstring AIX,$(BFARCH)))
else
ifneq (,$(findstring MacOSX,$(BFARCH)))
   override DEFOPT += -noObjy
else
ifneq (,$(findstring HP-UX,$(BFARCH)))
##   override DEFOPT += -noObjy-aCC
	# Overridden by gg 25/03/98
##   override DEFOPT += -Objy-aCC-noOptimize
	## Tuned by gg 14/12/98
	override DEFOPT := $(subst native,aCC,$(DEFOPT))
	ifeq (,$(findstring aCC,$(DEFOPT)))
		DEFOPT += -aCC
	endif
else
ifneq (,$(findstring IRIX,$(BFARCH)))
else
ifneq (,$(findstring OSFV,$(BFARCH)))
else
ifneq (,$(findstring SunOS,$(BFARCH)))

endif   # AIX
endif   # MacOSX
endif   # HP-UX
endif   # IRIX
endif   # OSFV
endif   # SunOS

#__________________________________________________________________________
#++ arrange switches for shared and static libraries.

    # set -Lstatic as default
    override DEFOPT += -Lstatic

    #-> correct common typo such as -LShared
    ifneq (,$(findstring -LShared,$(ROPT)))
	override ROPT += -noLShared-Lshared
    endif
    ifneq (,$(findstring -LStatic,$(ROPT)))
	override ROPT += -noLStatic-Lstatic
    endif

    #-> Static and Shared are exclusive	
    ifneq (,$(findstring -Static,$(ROPT)))
	ifneq (,$(findstring -Shared,$(ROPT)))
	    ERROR:=$(shell echo Error: use either -Shared or -Static >& 2)
	    abort
	endif
    endif

    # set default for ROPT= or ROPT=-Static or ROPT=-Shared
    ifeq (,$(findstring Lstatic,$(ROPT)))
	ifeq (,$(findstring Lshared,$(ROPT)))
	    ifeq (,$(findstring Static,$(ROPT)))
		ifneq (,$(findstring Shared,$(ROPT)))
		    override ROPT += -Shared-Lstatic-Lshared
		else
		    override ROPT += -Static-Lstatic
		endif
	    else
		ifeq (,$(findstring -Shared,$(ROPT)))
		    override ROPT += -Static-Lstatic
		endif
	    endif	
	endif
    endif

#__________________________________________________________________________
#++ setting PROPT, default of production release

ifneq (,$(findstring Production,$(ROPT)))
    export PRODUCTION_RELEASE := 1
endif

ifneq (,$(PRODUCTION_RELEASE))
   override PROPT := -noDebug-Prunebin-Time-Fastbuild
ifneq (,$(findstring AIX,$(BFARCH)))
else
ifneq (,$(findstring HP-UX,$(BFARCH)))
else
ifneq (,$(findstring IRIX,$(BFARCH)))
else
ifneq (,$(findstring OSFV,$(BFARCH)))
else
ifneq (,$(findstring SunOS,$(BFARCH)))
endif
endif
endif
endif
endif
endif  # PRODUCTION_RELEASE

#__________________________________________________________________________
#++ OPTIONS = DEFOPT + BFARCH + PROPT + ROPT [1,2,3]
#   [1]. DEFOPT is default for all users.
#        PROPT is default for production release.
#	 ROPT is command-line options, normally for Release Coordinator. 
#   [2]. Put double word architecture (e.g. HP-UX) into 'longarch'.
#   [3]. This section will keep only on-switch in OPTIONS.
export DEFOPT PROPT OPTIONS
longarch := HP ADDMORE
override arch := $(shell perl -e '$$arch=$$ARGV[0]; shift; $$longarch=join(" ",@ARGV); split(/\-/, $$arch); if ($$longarch=~/.*$$_[0].*/) {print "$$_[0]-$$_[1]"} else {print "$$_[0]"}' $(BFARCH) $(longarch))
override suffix := $(subst $(arch),,$(BFARCH))
override OPTIONS :="$(DEFOPT) $(PROPT) $(suffix) $(ROPT)"

ifdef ROPT
    ifneq (-,$(shell echo $(shell echo $(ROPT) | cut -c1)))
	ERROR := $(shell echo "error: First letter in ROPT should be '-'" >&2) 
	exit 1
    endif
endif

#__________________________________________________________________________
#++ function to clean up the OPTIONS. It will be called twice, the first
#   pass filter out all the no<option> switches. The second pass put 
#   switchs in proper order.

define CLEANOPT
perl -e '$$options = join (" ", @ARGV); $$options=~s/ //g; split(/\-/, $$options); shift @_; $$i=0; foreach $$opt (@_) {($$var)=$$opt=~/no(.*)/; if ("$$var" eq "") {$$list[$$i]=$$opt; $$listq[$$i++]=1;} else {for ($$j=0; $$j <= $$#list; $$j++) {if ($$list[$$j] eq "$$var"){$$listq[$$j]=0;}}}}; for ($$i=0; $$i <= $$#list; $$i++) {eval "\$$def=\$$$$list[$$i]";if ($$listq[$$i] == 1 && $$def == 0) {print "-$$list[$$i]";eval "\$$$$list[$$i]=1"}}' pad $(OPTIONS)
endef

override OPTIONS := $(shell $(CLEANOPT))

ifdef DEBUG
ERROR := $(shell echo "[1] DEFOPT=$(DEFOPT)\nPROPT=$(PROPT)\nROPT=$(ROPT)\narch=$(arch)\nOPTIONS-pass1=$(OPTIONS)" >& 2)
endif

#__________________________________________________________________________
#++ turn on -Skipdfile for -Skiplib
ifneq (,$(findstring -Skiplib,$(OPTIONS)))
    override OPTIONS += -Skipdfile
endif	

#__________________________________________________________________________
#++ temporary implementation
ifneq (,$(findstring -NOBIN,$(OPTIONS)))
    export NOBIN := 1
endif	
ifdef NOBIN
    override OPTIONS += -NOBIN
endif

#__________________________________________________________________________
#++ -gcc and -native are exclusive 
ifneq (,$(findstring -gcc,$(OPTIONS)))
override OPTIONS += -nonative
endif

#__________________________________________________________________________
#-> temporary switch -cxx5 for CXX5.0
ifneq (,$(findstring OSF1V4,$(BFARCH)))
    ifneq (,$(shell cxx -V 2>&1 | grep "C++ V5"))
	override OPTIONS += -cxx5
    endif
endif

#__________________________________________________________________________
#++ define the proper order of switchs in OPTIONS
#   first option comes first.	

    ifneq (,$(findstring -Static,$(OPTIONS)))
	override OPTIONS += -noStatic-Static
    endif
    ifneq (,$(findstring -Shared,$(OPTIONS)))
	override OPTIONS += -noShared-Shared
    endif
    ifneq (,$(findstring -Lshared,$(OPTIONS)))
	override OPTIONS += -noLshared-Lshared
    endif
    ifneq (,$(findstring -Lstatic,$(OPTIONS)))
	override OPTIONS += -noLstatic-Lstatic
    endif

ifdef DEBUG
ERROR := $(shell echo "\n[2] OPTIONS=$(OPTIONS)\n" >& 2)
endif

override OPTIONS := $(shell $(CLEANOPT))
override OPTIONS :="$(arch)$(OPTIONS)"

#__________________________________________________________________________
# Set ROOTSYS if it is not defined by srtpath

ifndef ROOTSYS
ifeq ($(MAKECMDGOALS), bbROOT)
    PKG = SoftRelTools
    MK=arch_spec_ROOT.mk
    include $(shell $(findmk))

bbROOT:
	echo "$(BBRROOTSYS)"

else
    export ROOTSYS := $(shell gmake bbROOT NOTFIRST=1)
endif   # bbROOT
endif   # no ROOTSYS

endif  # FIRST TIME

#___________________________________________________________________________
#++ set -Ldlink for ldlink, bootstrap, siteinstall

#ifneq (,$(findstring ldlink,$(MAKECMDGOALS)))
#    override OPTIONS+=-Ldlink2
#endif

#___________________________________________________________________________
#++ list of database directory to be protected

export protect_dbdir := databases users objy Releases

#==========================================================================
#++ process the environment variables and gmake variables.

#__________________________________________________________________________
# Set up environment for ROPT=-Profile build.

ifneq (,$(findstring -Profile,$(OPTIONS)))
  override CXXFLAGS += -pg
  override CCFLAGS += -pg
  override FCFLAGS += -pg
  override FFCFLAGS += -pg
  override F90CFLAGS += -pg
  export CXXFLAGS CCFLAGS FCFLAGS FFCFLAGS F90CFLAGS
endif

#__________________________________________________________________________
# make the correct BFCURRENT value

newCURRENT := $(shell touch .current; cat .current)
ifneq ($(newCURRENT),)
  # if the .current file exists in the release directory, use it
  ifneq ($(BFCURRENT),$(newCURRENT))
     ERROR := $(shell echo This release based on $(shell cat .current), overriding BFCURRENT value of $(BFCURRENT) >& 2)
  endif
  override BFCURRENT := $(newCURRENT)
  export BFCURRENT
endif

#__________________________________________________________________________
#++ include list of special packages for BaBar

PKG=SoftRelTools
MK=babar_specials.mk
export INMAIN := 1 
include $(shell $(findmk))
INMAIN :=

#__________________________________________________________________________
# set up search order - first the tmp local tree, tmp current, etc

tmparch :=
topdir :=

ifdef BFCURRENT  
ifneq (,$(findstring -Shared,$(OPTIONS)))
tmparch :=-I$(BFDIST)/releases/$(BFCURRENT)/shtmp/$(BFARCH) -I$(BFDIST)/releases/$(BFCURRENT)/tmp/$(BFARCH)
else
tmparch :=-I$(BFDIST)/releases/$(BFCURRENT)/tmp/$(BFARCH)
endif
topdir :=-I$(BFDIST)/releases/$(BFCURRENT)

# add customized system header files to search path
ifneq (,$(findstring HP-UX,$(BFARCH)))
Bfarch=$(word 1,$(subst -, ,$(BFARCH)))-$(word 2,$(subst -, ,$(BFARCH)))
else
Bfarch=$(word 1,$(subst -, ,$(BFARCH)))
endif
endif

ifneq (,$(findstring -Shared,$(OPTIONS)))
override CPPFLAGS += -I$(TOPDIR)/shtmp/$(BFARCH) $(tmparch)
else
override CPPFLAGS += -I$(TOPDIR)/tmp/$(BFARCH) $(tmparch)
endif
override CPPFLAGS += -I$(TOPDIR) $(topdir)

# Include file lookup path when processing IDL files is the cwd, then
# the top of the test release, then the top of the base release
override IDLINCLFLAGS += -I$(TOPDIR) $(topdir)

ifndef BFOVERRIDE
BFOVERRIDE = .
endif

ifneq (,$(findstring Lshared,$(OPTIONS)))
override LDFLAGS +=  -L$(TOPDIR)/shlib/$(BFARCH) -L$(BFDIST)/releases/$(BFCURRENT)/shlib/$(BFARCH)
endif

# the following four assignments are overrides because they matter so much,
# yet users might want to include various options
override CPPFLAGS += -I$(TOPDIR)/include
ifneq (,$(findstring Lstatic,$(OPTIONS)))
override LDFLAGS += -L$(TOPDIR)/lib/$(BFARCH)
endif

ifdef BFCURRENT  
override CPPFLAGS += -I$(BFDIST)/releases/$(BFCURRENT)/include
ifneq (,$(findstring Lstatic,$(OPTIONS)))
override LDFLAGS += -L$(BFDIST)/releases/$(BFCURRENT)/lib/$(BFARCH)
endif
endif

override libdir = $(TOPDIR)/lib/$(BFARCH)/
override shlibdir = $(TOPDIR)/shlib/$(BFARCH)/
override bindir = $(TOPDIR)/bin/$(BFARCH)/
override includedir = $(TOPDIR)/include/
override workdir=$(TOPDIR)/tmp/$(BFARCH)/
override shworkdir=$(TOPDIR)/shtmp/$(BFARCH)/
override mandir = $(TOPDIR)/man/
override docdir = $(TOPDIR)/doc/
override testdir = $(TOPDIR)/test/$(BFARCH)/
override dbdir = $(TOPDIR)/database

ifneq (,$(findstring Shared,$(OPTIONS)))
   export templatedir := $(shlibdir)templates/
else
   export templatedir := $(libdir)templates/
endif

# force local executables onto end of PATH
export PATH := $(TOPDIR)/bin/$(BFARCH):$(PATH)

#++ define schema build sequence here. The sequence is basically determined
# by the PackageList link sequence in reverse order.

PKG = PackageList
MK=link_all_reco.mk
include $(shell $(findmk))
ifneq (,$(DUMP_PACKAGELIST))
ERROR := $(shell echo PACKAGELIST: $(PACKAGELIST) >& 2)
endif

#___________________________________________________________________________
# You might want to override the following lines
# with an explicit list in your own, private makefiles

ifndef PACKAGES
ifdef DEBUG11
  ERROR := $(shell echo "define variable PACKAGES [DEBUG11-1]" >& 2)
endif
ifneq (1,$(shell if [ -d .online ]; then echo 1; fi))
    export PACKAGES := $(subst /,, $(dir $(wildcard */GNUmakefile)))
else   # online
    export PACKAGES := $(filter-out \
		$(subst /,, $(dir $(wildcard */Makefile))), \
		$(subst /,, $(dir $(wildcard */GNUmakefile))))
endif  # online
    ifeq (,$(findstring database,$(PACKAGES)))
         override PACKAGES += database   # add database to PACKAGES
    endif
endif    # PACKAGES

#___________________________________________________________________________
#++ exclude package from list if it is on only_sunos5 list.

ifeq (,$(findstring SunOS5,$(BFARCH)))
    export PACKAGES := $(filter-out $(only_sunos5), $(PACKAGES))
endif

#___________________________________________________________________________
# The following is a list of packages in the release that have DDL (.ddl)
# files.  Note that the objy_specials packages don't conform to the BaBar
# standard for placement of .ddl files and therefore must be added
# explicitly if they exist.

ifndef SCHEMAPACK
export SCHEMAPACK := $(subst /,, $(sort $(foreach var, $(PACKAGES), $(dir $(wildcard $(var)/*.?dl)))))
SCHEMAPACK += $(filter $(PACKAGES),$(objy_specials))
endif    # SCHEMAPACK

ifndef SCHEMALIST
export SCHEMALIST  = $(filter $(SCHEMAPACK),$(PACKAGELIST))
endif    # SCHEMALIST

#==========================================================================
#+++ Targets

#-> equivalent to 'lib', followed by 'bin', this is the default target
ifeq (,$(findstring Linux2,$(BFARCH)))
	SPECNUMLOCATION := 3
else
	SPECNUMLOCATION := 4
endif

#___________________________________________________________________________
# The 'all' target. It is called by the 'cold' target.

all:
	-@$(slaclog1)
	@echo "-> all:"
	-@newest=`ls -t .spec.$(BFARCH).* 2>/dev/null | head -1`; \
	 if [ "$$newest" ]; then OIFS=$$IFS; IFS=.; set $$newest; \
	    IFS=$$OIFS; spec=".spec.$(BFARCH).`expr $$$(SPECNUMLOCATION) + 1`"; \
	    else spec=".spec.$(BFARCH).1"; fi; \
	 echo "Saving Spec of this build to $$spec"; \
	 $(MAKE) spec 1>$$spec; \
	 cp $$spec .spec-$(BFARCH)
	@echo "Start saving index of binaries into bin/$(BFARCH)/Index"
	@$(RM) -f bin/$(BFARCH)/Index	
	@echo "Index of $(BFARCH) binaries:" > bin/$(BFARCH)/Index
	@$(MAKE) installdirs
ifeq (,$(PRODUCTION_RELEASE)$(PRODUCTION_TEST))
  ifeq (,$(BYPASS_DB_IMPORT))
	@$(MAKE) database.import;
  endif
endif
ifeq (,$(BYPASS_SCHEMA))
	@$(MAKE) schema
endif

ifdef IS_COLD_BUILD
ifneq (,$(findstring SunOS,$(BFARCH)))
	@$(MAKE) fastinclude
else
ifneq (,$(findstring Linux2,$(BFARCH)))
	@$(MAKE) fastinclude
else
	@${MAKE} include
endif
endif
else
	@${MAKE} include
endif

	@$(MAKE) debugger

ifdef IS_COLD_BUILD
ifneq (,$(findstring SunOS,$(BFARCH)))
	@$(MAKE) fastlib
else
ifneq (,$(findstring Linux2,$(BFARCH)))
	@$(MAKE) fastlib
else
	@${MAKE} lib
endif
endif
else
	@${MAKE} lib
endif

ifneq (,$(PRODUCTION_RELEASE)$(PRODUCTION_TEST))
	-MULTIROOTPKGS="$(subst /,, $(dir $(wildcard */root_*)))"; \
	if [ -n "$$MULTIROOTPKGS" ]; then \
	    echo "============ Creating multiple package ROOT shared libraries: Start ============="; \
	    for pkg in $$MULTIROOTPKGS; do \
		$(MAKE) $$pkg.root; \
	    done; \
	    echo "============ Creating multiple package ROOT shared libraries: End ============="; \
	fi
endif  # PRODUCTION_RELEASE

ifdef NOBIN
ifdef IS_COLD_BUILD
ifneq (,$(findstring SunOS,$(BFARCH)))
	@$(MAKE) fastbincomp IS_COLD_BUILD=yes
else
ifneq (,$(findstring Linux2,$(BFARCH)))
	@$(MAKE) fastbincomp IS_COLD_BUILD=yes
else
	@${MAKE} bincomp IS_COLD_BUILD=yes
endif
endif
else
	@${MAKE} bincomp
endif
endif

ifneq (,$(PRODUCTION_RELEASE)$(PRODUCTION_TEST))
	@echo "============ Creating and Testing Database Loaders: Start ============="
	@if test -d BdbEvsLoaders;     then $(MAKE) BdbEvsLoaders.bin;    fi;
	@if test -d BdbCondLoaders;    then $(MAKE) BdbCondLoaders.bin;   fi;
ifneq (,$(TEST_DATABASE_RESET))
	@$(MAKE) database.reset
endif
	@echo "============ Creating and Testing Database Loaders: End ============="
endif

ifneq (,$(PRODUCTION_RELEASE)$(PRODUCTION_TEST))
	@echo "============ Linking special executables for the nightly: Start ============="
	@if test -d BdbAdminTools;     then $(MAKE) BdbAdminTools.BdbDomainBootNames; fi
	@if test -d CdbTools;	       then $(MAKE) NOBIN= CdbTools.bin; fi
	@if test -d Kanga;	       then $(MAKE) Kanga.Kanga; fi
	@echo "============ Linking special executables for the nightly: End ============="
endif

ifndef NOBIN
ifdef IS_COLD_BUILD
ifneq (,$(findstring SunOS,$(BFARCH)))
	@$(MAKE) fastbin IS_COLD_BUILD=yes
else
ifneq (,$(findstring Linux2,$(BFARCH)))
	@$(MAKE) fastbin IS_COLD_BUILD=yes
else
	@${MAKE} bin IS_COLD_BUILD=yes
endif
endif
else
	@${MAKE} bin
endif
endif

ifneq (,$(findstring Prunebin,$(OPTIONS)))
	cd $(bindir); \
	for f in *; do \
	    if [ ! -s $$f ] ;then \
		rm $$f; \
	    fi; \
	done
endif
	-@$(MAKE) java
ifneq (,$(findstring Linux24,$(BFARCH)))
  ifeq (,$(findstring Optimize,$(BFARCH)))
    ifneq (,$(PRODUCTION_RELEASE)$(PRODUCTION_TEST))
	@$(MAKE) glimpse
	-@echo "find . -type d -exec chgrp bfadmin {} \\;" ; \
	find . -type d -exec chgrp bfadmin {} \;
	-@if [ -d "tmp/$(BFARCH)/java" ]; then \
	    echo "find tmp/$(BFARCH)/java -exec chmod a+r {} \\;" ; \
	    find tmp/$(BFARCH)/java -exec chmod a+r {} \; ; \
	fi
	-@if test -d Kanga; then $(MAKE) Kanga.kangaDoc; fi;
    endif
  endif
endif
	-@$(slaclog2)

#__________________________________________________________________________
#++ bootstrap of SRT, build package in specific order
#   Note: bootstrap builds two types of targets:
#     1. a target that builds critical include/library/binary 
#     2. a target that must be built first before other packages. Add prefix '-' so that
#	 it won't abort the production release.  

bootstrap:
	-@$(slaclog1)
	@echo "BOOTSTRAP: build basic packages in specific order"
	@$(MAKE) installdirs
	@$(MAKE) $(OVERRIDES) ldlink
ifneq (,$(PRODUCTION_RELEASE)$(PRODUCTION_TEST))
  ifneq (,$(REQUIRE_DB_EXISTS))
	@$(MAKE) database.exists
	echo "Removing database DDL generated files for this architecture."
	$(RM) -f $(TOPDIR)/tmp/$(BFARCH)/*/*_ddl.*
  else
	@$(MAKE) database.cleanup
    ifeq (,$(REQUIRE_DB_IMPORT))
      ifeq (,$(BYPASS_DB_DELETEBOOT)$(REBUILD))
	@$(MAKE) database.deleteboot CONFIRM_DELETEBOOT=yes
      endif
      ifeq (,$(BYPASS_DB_MAKEBOOT)$(REBUILD))
        ifeq (,$(BYPASS_REFERENCE_FDB))
	@$(MAKE) database.import IMPORT_REFERENCE_FDB=yes REQUIRE_DB_IMPORT=yes BYPASS_LOAD=yes NFS_CACHE_WAIT=yes
        else
	@$(MAKE) database.makeboot CONFIRM_MAKEBOOT=yes NFS_CACHE_WAIT=yes
        endif
      endif
    else
      ifeq (,$(BYPASS_DB_DELETEBOOT))
	@$(MAKE) database.deleteboot CONFIRM_DELETEBOOT=yes
      endif
	@$(MAKE) database.import BYPASS_LOAD=yes NFS_CACHE_WAIT=yes
    endif
  endif
endif
	@if test -d SoftRelTools; then $(MAKE) SoftRelTools.all; fi;
	@if test -d geant4; then $(MAKE) geant4.include; fi;
	-if test -d RegrTest; then $(MAKE) RegrTest.all; fi;
	-if test -d beget; then $(MAKE) beget.lib; fi;
	-if test -d OrcSmi; then $(MAKE) OrcSmi.binscripts; fi;
	@echo "_______END of BOOTSTRAP________________________________________"
	-@$(slaclog2)

#__________________________________________________________________________
#++ target cold is used when making a production release
cold:
	-@$(slaclog1) 
	@echo "-> cold:"
	@echo "_______________________________________________________________"
	@echo "      Starting COLD build of BaBar Software Release $(BFCURRENT)"
	@echo "                       "DATE = `date` 
	@echo "                       "HOST = `hostname`
	@echo "                       ""Operating System = $(BFARCH)"
	@echo "                       ""Package Version = $(BFCURRENT)"
	@echo "                       ""BFOVERRIDE = $(BFOVERRIDE)"
	@echo "                       ""BFDIST     = $(BFDIST)"
	@echo "                       ""HPBUILD    = $(HPBUILD)"
	@echo "                       ""AIXBUILD   = $(AIXBUILD)"
	@echo "                       ""OSFBUILD   = $(OSFBUILD)"
	@echo "                       ""BFSITE     = $(BFSITE)"
	@echo "_______________________________________________________________"
	@$(MAKE) NOBIN= bootstrap
	@$(MAKE) all IS_COLD_BUILD=yes
#	@if test -d gnbase; then $(MAKE) gnbase.bin MODE=INTERACTIVE; fi;
	-if test -d XrdSeckrb4; then $(MAKE) XrdSeckrb4.libext; fi;
	-if test -d XrdSeckrb5; then $(MAKE) XrdSeckrb5.libext; fi;
	-if test -d XrdSecgsi; then $(MAKE) XrdSecgsi.libext; fi;
	-if test -d XrdOfs; then $(MAKE) XrdOfs.libext; fi;
	-if test -d XrdRootd; then $(MAKE) XrdRootd.libext; fi;
	-if test -d RooFitCore -a -r RooFitCore/root_RooFitCore.mk; then $(MAKE) RooFitCore.root ROPT+=-Shared; fi;
	-if test -d RooFitModels -a -r RooFitModels/root_RooFitModels.mk; then $(MAKE) RooFitModels.root ROPT+=-Shared; fi;
	-if test -d RooFitBabar -a -r RooFitBabar/root_RooFitBabar.mk; then $(MAKE) RooFitBabar.root ROPT+=-Shared; fi;
	-@$(slaclog2)

#__________________________________________________________________________
#-> nightly built, cold target minus bin
nightly:
	echo "-> nightly:"
	@$(MAKE) cold NOBIN=1 NOTFIRST=

#__________________________________________________________________________
#-> Full production release
production:
	@$(MAKE) cold PRODUCTION_RELEASE=yes NOTFIRST=

# the following two rules are needed because the existing directories
# confuse the more general rule

.PHONY: schema include lib bin test debugger .dbxinit .gdbinit

debugger: .dbxinit .gdbinit

#__________________________________________________________________________
.dbxinit:
	-if [ -w ./ ] ; then \
	packages="$(PACKAGES)"; \
	if [ -f $@ ] ; then \
	verb=Updating; else verb=Creating; fi; \
	echo $$verb $< debugger directory file $@; \
	echo "use \\" > $@; \
	for v in $$packages; do \
	v="$(TOPDIR)/$$v"; \
	for vv in `find $$v -type d`; do \
	if [ `basename $$vv` != "CVS" ]; then echo " $$vv \\" >> $@; fi;\
	done; done; \
	for v in $$packages; do \
	v="$(TOPDIR)/tmp/$(BFARCH)/$$v"; \
	if [ -d $$v ]; then echo " $$v \\" >> $@ ; fi; done; \
	echo " $(TOPDIR)/include" >> $@; \
	if [ -w gnbase/ ] ; then \
	cd gnbase; $(RM) $@ > /dev/null 2>&1; \
	ln -s ../$@ . > /dev/null 2>&1; fi; fi

#__________________________________________________________________________
.gdbinit:
	-if [ -w ./ ] ; then \
	packages="$(PACKAGES)"; \
	if [ -f $@ ] ; then \
	verb=Updating; else verb=Creating; fi; \
	echo $$verb $< debugger directory file $@; \
	$(RM) $@ > /dev/null 2>&1; \
	for v in $$packages; do \
	v="$(TOPDIR)/$$v"; \
	for vv in `find $$v -type d`; do \
	if [ `basename $$vv` != "CVS" ]; then echo dir $$vv >> $@; fi; \
	done; done; \
	for v in $$packages; do \
	v="$(TOPDIR)/tmp/$(BFARCH)/$$v"; \
	if [ -d $$v ]; then echo dir $$v >> $@ ; fi; done; \
	echo "dir $(TOPDIR)" >> $@; \
	if [ -w gnbase/ ] ; then cd gnbase; \
	$(RM) $@ > /dev/null 2>&1; \
	ln -s ../$@ . > /dev/null 2>&1; fi; fi

#__________________________________________________________________________
#++ bin, lib, etc

# This is a macro to aid in doing stuff to all packages
define buildAll
	for package in xx $(filter-out $(BUILDALLSKIP_$@),$(PACKAGES)); do \
	    if [ "$$package" != "xx" -a "$$package" != "BASE" ]; then \
		$(MAKE) $(BUILDALLMAKEFLAGS_$@) $$package.$@ NOPKG=yes; \
	    fi; \
	done
endef

BUILDALLSKIP_bin := $(skip_bin)
BUILDALLMAKEFLAGS_bin := LOOPPKG=1
bin:
	-@$(slaclog1) 
	@echo "-> bin"; 
	$(buildAll)
ifdef IS_COLD_BUILD
	@echo "-> rootFileMap";
	mapRootFiles -c 1> tmp/$(BFARCH)/SoftRelTools/rootFileMap
endif
	-@$(slaclog2) 

ifeq (,$(findstring OSF,$(BFARCH)))
export BINPACKAGES := $(filter-out $(BUILDALLSKIP_bin),$(PACKAGES))
export BINPACKAGES_1 := $(wordlist 1,300,$(BINPACKAGES))
export BINPACKAGES_2 := $(wordlist 301,600,$(BINPACKAGES))
export BINPACKAGES_3 := $(wordlist 601,900,$(BINPACKAGES))
export BINPACKAGES_4 := $(wordlist 901,$(words $(BINPACKAGES)),$(BINPACKAGES))
endif
export BINOPTIONS := NOPKG=yes
export MAKE
fastbin:
	-@$(slaclog1) 
	@echo "-> fastbin";
	makeParallelBin;
ifdef IS_COLD_BUILD
	@echo "-> rootFileMap";
	mapRootFiles -c 1> tmp/$(BFARCH)/SoftRelTools/rootFileMap
endif
	-@$(slaclog2) 

ifeq (,$(findstring OSF,$(BFARCH)))
export BINCOMPPACKAGES := $(filter-out $(BUILDALLSKIP_bincomp),$(PACKAGES))
export BINCOMPPACKAGES_1 := $(wordlist 1,300,$(BINCOMPPACKAGES))
export BINCOMPPACKAGES_2 := $(wordlist 301,600,$(BINCOMPPACKAGES))
export BINCOMPPACKAGES_3 := $(wordlist 601,900,$(BINCOMPPACKAGES))
export BINCOMPPACKAGES_4 := $(wordlist 901,$(words $(BINCOMPPACKAGES)),$(BINCOMPPACKAGES))
endif
export BINCOMPOPTIONS := NOPKG=yes
export MAKE
fastbincomp:
	-@$(slaclog1) 
	@echo "-> fastbincomp";
	makeParallelBincomp;
ifdef IS_COLD_BUILD
	@echo "-> rootFileMap";
	mapRootFiles -c 1> tmp/$(BFARCH)/SoftRelTools/rootFileMap
endif
	-@$(slaclog2) 

BUILDALLSKIP_testbin := $(skip_bin)
testbin:
	-@$(slaclog1) 
	@echo "-> testbin";
	$(buildAll)
	-@$(slaclog2) 

BUILDALLSKIP_allbin := $(skip_bin)
allbin:
	-@$(slaclog1) 
	@echo "-> allbin";
	$(buildAll)
	-@$(slaclog2) 

BUILDALLSKIP_extrabin := $(skip_bin)
extrabin:
	-@$(slaclog1) 
	@echo "-> extrabin";
	$(buildAll)
	-@$(slaclog2) 

BUILDALLSKIP_binscripts := $(skip_bin)
binscripts:
	-@$(slaclog1) 
	@echo "-> binscripts";
	$(buildAll)
	-@$(slaclog2) 

BUILDALLSKIP_extrabinscripts := $(skip_bin)
extrabinscripts:
	-@$(slaclog1) 
	@echo "-> extrabinscripts";
	$(buildAll)
	-@$(slaclog2) 

BUILDALLSKIP_lib := $(skip_lib)
lib:
	-@$(slaclog1) 
	@echo "-> lib";
	$(buildAll)
	-@$(slaclog2) 

BUILDALLSKIP_libext := $(skip_lib)
libext:
	-@$(slaclog1) 
	@echo "-> libext";
	$(buildAll)
	-@$(slaclog2) 

BUILDALLSKIP_bincomp := $(skip_bin)
bincomp:
	-@$(slaclog1) 
	@echo "-> bincomp";
	$(buildAll)
ifdef IS_COLD_BUILD
	@echo "-> rootFileMap";
	mapRootFiles -c 1> tmp/$(BFARCH)/SoftRelTools/rootFileMap
endif
	-@$(slaclog2) 

ifeq (,$(findstring OSF,$(BFARCH)))
LIBPACKAGES := $(filter-out $(BUILDALLSKIP_lib),$(PACKAGES))
export LIBPACKAGES_1 := $(wordlist 1,300,$(LIBPACKAGES))
export LIBPACKAGES_2 := $(wordlist 301,600,$(LIBPACKAGES))
export LIBPACKAGES_3 := $(wordlist 601,900,$(LIBPACKAGES))
export LIBPACKAGES_4 := $(wordlist 901,$(words $(LIBPACKAGES)),$(LIBPACKAGES))
endif
export LIBOPTIONS := NOPKG=yes
export MAKE
fastlib:
	-@$(slaclog1) 
	@echo "-> fastlib";
	makeParallelLib;
	-@$(slaclog2) 

BUILDALLSKIP_test := $(skip_test)
test:
	-@$(slaclog1) 
	@echo "-> test";
	$(buildAll)
	-@$(slaclog2)

siteinstall:
	if [ -w SoftRelTools ]; then \
		$(MAKE) SoftRelTools.binscripts SITEINSTALL=1; \
	fi
	$(MAKE) binscripts SITEINSTALL=1
	$(MAKE) python SITEINSTALL=1
	$(MAKE) ldlink SITEINSTALL=1

#######################################################################
#++ ldlink: make ld links for shared SYSLIB

ifeq ($(MAKECMDGOALS), ldlink)    # ldlink
    export LDLINKSUFFIX := _LDLINK

#_______________________________________________________________________
#++ Include ldlink.mk

PKG=SoftRelTools
MK=ldlink.mk
include $(shell $(findmk))

#__________________________________________________________________________
#++ Include SoftRelTools/arch_spec_*.mk and their add-on mk-, mk+ files.

ifndef DEBUG28
    MKFILES := $(notdir $(wildcard $(SRTPATH)/arch_spec_*.mk) arch_spec.mk)
    MKFILES2 := $(MKFILES:.mk=.mk-)
    MKFILES3 := $(MKFILES:.mk=.mk+)
    -include $(foreach MK, $(MKFILES2), $(shell $(findmk))) $(foreach MK, $(MKFILES), $(shell $(findmk))) $(foreach MK, $(MKFILES3), $(shell $(findmk)))
endif   # DEBUG28

#__________________________________________________________________________
#++ ldlink: make ld links for shared SYSLIB

BUILDALLSKIP_lib :=
ldlink:
	-@$(slaclog1)
	@echo "-> ldlink:";
ifndef DEBUG28
	$(ldlink)
else
	$(buildAll)
endif  # DEBUG28
	# Need to make a symlink for some Intel compiler libraries for runtime
ifneq (,$(findstring icc,$(BFARCH)))
	for file in $(ICC_DIR)/lib/lib*.so* $(IFORT_DIR)/lib/lib*.so*; do \
		$(RM) -f shlib/$(BFARCH)/ldlink/`basename $$file`; \
		ln -s $$file shlib/$(BFARCH)/ldlink; \
		icount=`expr $$icount + 1`; \
	done; \
	echo "make $$icount links for Intel Compiler"
endif # Intel CC
	-@$(slaclog2)

endif # ldlink

#######################################################################

#__________________________________________________________________________
#++ define doxygen

PKG=SoftRelTools
MK=arch_spec_Doxygen.mk
-include $(shell $(findmk))
PKG=SoftRelTools
MK=arch_spec_PERL.mk
-include $(shell $(findmk))

doxygen:
	-@$(slaclog1)
	@echo "-> doxygen"
ifeq (Linux2,$(BFARCH))
	sed 's%@PERL@%$(PERL)%g; \
		s:@RELEASE@:$(BFDIST)/releases/$(BFCURRENT):g; \
		s:@BFCURRENT@:$(BFCURRENT):g' SoftRelTools/Doxyfile >Doxyfile
	$(DOXYGEN)
	rm -f Doxyfile
else
	@echo "Only done on Linux2"
endif
	-@$(slaclog2)

#__________________________________________________________________________
#++ interactive target.
interactive:
	-@$(slaclog1)
	if [ x$(findstring Inter,$(BFARCH))x = xx ]; then \
		echo "Not building interactive build (should have Inter in BFARCH)" ;\
		exit 1; \
	fi
	if [ x$(findstring Shared,$(ROPT))x = xx ]; then \
		echo "Not building with shared libraries (should include Shared in ROPT)" ; \
		exit 1 ; \
	fi
	$(MAKE) bootstrap LINK_ROOT=yes
	$(MAKE) include LINK_ROOT=yes
	-@if test -d RooUtils; then $(MAKE) RooUtils.binscripts LINK_rw=yes LINK_Tcl=yes; fi;
	-@if test -d BaBar; then $(MAKE) BaBar.lib LINK_ROOT=yes; fi;
	-@if test -d CLHEP; then $(MAKE) CLHEP.lib LINK_ROOT=yes; fi;
	-@if test -d BField; then $(MAKE) BField.lib LINK_ROOT=yes; fi;
	-@if test -d BtaEnv; then $(MAKE) BtaEnv.lib LINK_ROOT=yes; fi;
	-@if test -d ProxyDict; then $(MAKE) ProxyDict.lib LINK_ROOT=yes; fi;
	-@if test -d ErrLogger; then $(MAKE) ErrLogger.lib LINK_ROOT=yes; fi;
	-@if test -d AbsEvent; then $(MAKE) AbsEvent.lib LINK_ROOT=yes; fi;
	-@if test -d AbsCalo; then $(MAKE) AbsCalo.lib LINK_ROOT=yes; fi;
	-@if test -d AbsEnv; then $(MAKE) AbsEnv.lib LINK_ROOT=yes; fi;
	-@if test -d BbrGeom; then $(MAKE) BbrGeom.lib LINK_ROOT=yes; fi;
	-@if test -d CommonUtils; then $(MAKE) CommonUtils.lib LINK_ROOT=yes; fi;
	-@if test -d difAlgebra; then $(MAKE) difAlgebra.lib LINK_ROOT=yes; fi;
	-@if test -d TrajGeom; then $(MAKE) TrajGeom.lib LINK_ROOT=yes; fi;
	-@if test -d PepEnvData; then $(MAKE) PepEnvData.lib LINK_ROOT=yes; fi;
	-@if test -d TrkBase; then $(MAKE) TrkBase.lib LINK_ROOT=yes; fi;
	-@if test -d PDT; then $(MAKE) PDT.lib LINK_ROOT=yes; fi;
	-@if test -d ProbTools; then $(MAKE) ProbTools.lib LINK_ROOT=yes; fi;
	-@if test -d AbsPid; then $(MAKE) AbsPid.lib LINK_ROOT=yes; fi;
	-@if test -d BetaMicroAdapter; then $(MAKE) BetaMicroAdapter.lib LINK_ROOT=yes; fi;
	-@if test -d BetaRecoAdapter; then $(MAKE) BetaRecoAdapter.lib LINK_ROOT=yes; fi;
	-@if test -d Beta; then $(MAKE) Beta.lib LINK_ROOT=yes; fi;
	-@if test -d VtxBase; then $(MAKE) VtxBase.lib LINK_ROOT=yes; fi;
	-@if test -d BetaCoreTools; then $(MAKE) BetaCoreTools.lib LINK_ROOT=yes; fi;
	-@if test -d FrameUtil; then $(MAKE) FrameUtil.lib LINK_ROOT=yes; fi;
	-@if test -d Framework; then $(MAKE) Framework.lib LINK_ROOT=yes; fi;
	-@if test -d RooUtils; then $(MAKE) RooUtils.lib; fi;
	-@if test -d RooScribes; then $(MAKE) RooScribes.lib; fi;
	-@if test -d OdfCommon; then $(MAKE) OdfCommon.lib LINK_ROOT=yes; fi;
	-@if test -d BdbTime; then $(MAKE) BdbTime.lib LINK_ROOT=yes; fi;
	-@if test -d EidData; then $(MAKE) EidData.lib LINK_ROOT=yes; fi;
	-@if test -d EidDataR; then $(MAKE) EidDataR.lib; fi;
	-@if test -d TagData; then $(MAKE) TagData.lib LINK_ROOT=yes; fi;
	-@if test -d TagDataR; then $(MAKE) TagDataR.lib; fi;
	-@if test -d StdHepData; then $(MAKE) StdHepData.lib LINK_ROOT=yes; fi;
	-@if test -d StdHepDataR; then $(MAKE) StdHepDataR.lib; fi;
# This one is special b/c TrkContextEv adds a bunch of dependancies
	-@if test -d TrkFitter; then $(MAKE) TrkFitter.lib LINK_ROOT=yes LIBREMOVEFILES="TrkContextEv.cc TestPointHelix.cc"; fi;
	-@if test -d BtaCutTools; then $(MAKE) BtaCutTools.lib LINK_ROOT=yes; fi;
	-@if test -d BetaMicroBase; then $(MAKE) BetaMicroBase.lib LINK_ROOT=yes; fi;
	-@if test -d BetaMicroTruth; then $(MAKE) BetaMicroTruth.lib; fi;
	-@if test -d BtaCutImpl; then $(MAKE) BtaCutImpl.lib LINK_ROOT=yes; fi;
	-@if test -d BtaDataR; then $(MAKE) BtaDataR.lib; fi;
	-@if test -d GenEnv; then $(MAKE) GenEnv.lib LINK_ROOT=yes; fi;
	-@if test -d PepEnv; then $(MAKE) PepEnv.lib LINK_ROOT=yes; fi;
	-@if test -d RooCond; then $(MAKE) RooCond.lib; fi;
	-@if test -d RooModules; then $(MAKE) RooModules.lib; fi;
	-@if test -d BtaBuildEnv; then $(MAKE) BtaBuildEnv.lib LINK_ROOT=yes; fi;
	-@if test -d RooSequences; then $(MAKE) RooSequences.lib; fi;
	-@if test -d RooFramework; then $(MAKE) RooFramework.lib; fi;
	-@if test -d RooFitTools; then $(MAKE) RooFitTools.rootlib; fi;
	-@if test -d RooCutTools; then $(MAKE) RooCutTools.lib; fi;
	-@$(slaclog2) 

#__________________________________________________________________________
# root:

root:
	-@$(slaclog1) 
	@echo "-> bin"; 
	$(buildAll)
	-@$(slaclog2) 

#__________________________________________________________________________
# python:

BUILDALLSKIP_python := workdir database
python:
	-@$(slaclog1) 
	@echo "-> python"; 
	$(buildAll)
	-@$(slaclog2) 

#__________________________________________________________________________
#++ cleaning targets.

# clean is done by running cleanarch on all architectures

export BFARCHES :=  $(shell if [ "$(DEBUG15)" ]; then set -x; fi; \
	    list=; \
	    for v in xx $(wildcard bin/* lib/* shlib/* tmp/* shtmp/* test/*); do \
		if [ "$$v" != "xx" -a -d "$$v" ]; then \
			arch=`basename $$v`; \
		        if [ `expr " $$list " : ".* $$arch .*" ` -le 0 ]; then \
			    list="$$list `basename $$v`"; \
			fi; \
		fi; \
	    done; \
	    if [ -z "$$list" ]; then list=$(BFARCH); fi; \
	    echo $$list)

ALLLIST := bin/$(BFARCH) lib/$(BFARCH) shlib/$(BFARCH) $(shlibdir)ldlink $(shlibdir)LDLINK tmp/$(BFARCH) shtmp/$(BFARCH) test/$(BFARCH)
CLEANLIST := $(filter-out $(NOCLEAN:=%), $(ALLLIST))

clean:
	-@$(slaclog1) 
	@echo "-> clean:"
	@(true; $(foreach var,$(BFARCHES),$(MAKE) BFARCH=$(var) CLEANEXTRA=results cleanarch;))
	-@$(slaclog2) 

cleanbin:
	-@$(slaclog1) 
	@echo "-> cleanbin:"
	@echo $(RM) -r $(TOPDIR)/bin/$(BFARCH)
	-$(RM) -r $(TOPDIR)/bin/$(BFARCH)
	@echo mkdir $(TOPDIR)/bin/$(BFARCH)
	-mkdir $(TOPDIR)/bin/$(BFARCH)
	-@$(slaclog2) 

cleanlib:
	-@$(slaclog1) 
	@-ls -1da lib/$(BFARCH)/* 2>/dev/null | xargs $(RM) 2>/dev/null
	@echo "-> cleanlib:"
	-@$(slaclog2) 

# clean up incremental linking's object files
cleanincr:
	-@$(slaclog1) 
	@echo "-> cleanincr:"
	if [ "$(VERBOSE)" -eq "" ];then  echo $(RM) -r $(workdir)incr; fi
	$(RM) -r $(workdir)incr
	echo "remaking bin to make sure all binaries are up-to-date." 
	@$(MAKE) bin
	-@$(slaclog2) 

%.cleanincr:
	if [ "$(VERBOSE)" -eq "" ];then \
	  echo $(RM) -r $(workdir)incr/$(word 1,$(subst ., ,$@))/*.o; \
	fi
	$(RM) -r $(workdir)incr/$(word 1,$(subst ., ,$@))/*.o
	echo "remaking bin to make sure all binaries are up-to-date." 
	@$(MAKE) bin

# clean directories of one architecture.
cleanarch:
	-@$(slaclog1)
	@echo "-> cleanarch: $(BFARCH)"
	-list=; \
	if [ "$(DEBUG20)" ]; then set -x; fi; \
	for v in $(CLEANLIST) $(filter-out $(NOCLEAN:=%), $(CLEANEXTRA)) xx; do \
	    if [ "$vv" = "xx" ]; then continue; fi; \
	    dir="$$v"; \
	    if [ -h "$$v" ]; then \
	        linkto=`ls -ldt $$v | awk '{print $$NF}' -`; \
		if [ "$$linkto" = "." ]; then \
		   echo "'$$v' is soft link to '.', ignored"; \
		   dir=""; \
		elif [ -d "$$linkto" ]; then \
		   dir=`cd $$v; /bin/pwd`; \
		   echo "'$$v' is soft link to $$dir"; \
		fi; \
	    fi; \
	    if [ -d "$$dir" ]; then \
		list="$$list $$dir"; \
	    fi; \
	done; \
	if [ -n "$$list" ]; then \
	   echo "$(RM) -r $$list"; \
	   $(RM) -r $$list; \
	   echo "mkdir $$list"; \
	   mkdir $$list; \
	fi; 
	-$(MAKE) installdirs

	-@$(slaclog2)

# depend is done locally by explicit delete
# note this does not actually remake the .d files - that's
# done at the beginning of the next make
depend:
	@echo "-> depend:"
	$(RM) tmp/*/*/*.d

#__________________________________________________________________________
#++ JAVA variables and targets
#
JAVATARGETS = jar jclasses jclean jenv jdoc jbin jlib jcorba java javadoc javalinks
skip_include_d += $(JAVATARGETS)
.PHONY:	$(JAVATARGETS) 

skip_java += PARENT RELEASE database

jar jclasses jclean jenv jdoc jbin jlib :
	@for var in $(filter-out $(skip_java),$(PACKAGES)); do \
	   test -f $$var/java/GNUmakefile && $(MAKE) MAKEGOAL=$@ $$var.$@;\
	   done; true;

jcorba :
	packages=""; \
	for package in xx $(filter-out $(skip_java),$(PACKAGELIST)); do \
	    packages="$$package $$packages"; \
	done; \
	for package in $$packages; do \
	    if [ "$$package" != "xx" -a -f "$$package/java/GNUmakefile" ]; then \
		$(MAKE) MAKEGOAL=$@ $$package.$@ ; \
	    fi; \
	done;

java: javalinks jcorba jar jbin

PKG = SoftRelTools
MK=jstandard.mk
JAVAMK := $(shell $(findmk))

javadoc:
	@echo "-> javadoc:"
	$(MAKE) $(OVERRIDES) -f $(JAVAMK) fulldoc docdir=$(docdir)

javalinks:
	@echo "-> javalinks:"
	$(MAKE) $(OVERRIDES) -f $(JAVAMK) javalinks PACKAGES="$(PACKAGES)"

#__________________________________________________________________________
#++ Special doc variables and targets
#
DOCTARGETS = guide reference lightcpp glimpse
skip_include_d += $(DOCTARGETS)
.PHONY:	$(DOCTARGETS)

skip_guide += PARENT RELEASE 

guide:
	@(true; $(strip $(foreach var,$(filter-out $(skip_guide),$(PACKAGES)),test -f $(var)/doc/GNUmakefile && $(MAKE) MAKEGOAL=$@ $(var).$@; true;)))

PKG = SoftRelTools
MK=doctools.mk
DOCTOOLSMK:= $(shell $(findmk))

reference: 
	@echo "-> reference:"
	packages=""; \
	for package in xx $(filter-out $(skip_guide),$(PACKAGELIST)); do \
	    packages="$$package $$packages"; \
	done; \
	for package in $$packages; do \
	    if [ "$$package" != "xx" -a -d "$$package" ]; then \
		$(MAKE) -C $$package $(OVERRIDES) -f $(DOCTOOLSMK) reference PACKAGE=$$package ; \
	    fi; \
	done;

%.reference: 
	@echo "-> $@:   $(qtime)"; pkg=$(word 1,$(subst ., ,$@)); \
	$(MAKE) -C $(notdir $(basename $@)) $(OVERRIDES) -f $(DOCTOOLSMK) reference PACKAGE=$$pkg

glimpse: SoftRelTools.bin
	@echo "-> glimpse: building source index"
	srtglimpseindex > /dev/null 2>&1; \
	if [ ! -s .glimpse_partitions ]; then \
	    echo "warning: .glimpse_partitions is empty"; \
	fi		

#__________________________________________________________________________
#++ The schema target generates the database schema only. Note that the
#   packages must appear in a particular sequence determined by the inverted
#   PACKAGELIST list because of the DDL file processing.

schema: 
	-@$(slaclog1)
	@echo "-> schema:"
	packages=""; \
	for package in xx $(SCHEMALIST); do \
	    packages="$$package $$packages"; \
	done; \
	for package in $$packages; do \
	    if [ "$$package" != "xx" ]; then \
		$(MAKE) $$package.$@ NOPKG=yes; \
		if [ $$? -ne 0 -a "$(PRODUCTION_RELEASE)" != "" -a "$(SCHEMA_ERROR_IGNORE)" = "" ]; then \
		    exit 1; \
		fi; \
	    fi; \
	done;
ifneq (,$(PRODUCTION_RELEASE))
ifeq (,$(REQUIRE_DB_IMPORT)$(REQUIRE_DB_EXISTS)$(BYPASS_DB_EXPORT))
	$(MAKE) database.export FORCE_EXPORT=yes EXPORT_NAME=BaBarRef;
endif
endif
	-@$(slaclog2)

#__________________________________________________________________________
#++ The include target is routinely done for only those packages known to
#   require it. 

BUILDALLSKIP_include := $(skip_include) $(SCHEMAPACK)
include:
	@echo "-> include:"
	packages=""; \
	for package in xx $(SCHEMALIST); do \
	    packages="$$package $$packages"; \
	done; \
	for package in $$packages; do \
	    if [ "$$package" != "xx" ]; then \
		$(MAKE) $$package.$@ NOPKG=yes; \
		if [ $$? -ne 0 -a "$(PRODUCTION_RELEASE)" != "" -a "$(SCHEMA_ERROR_IGNORE)" = "" ]; then \
		    exit 1; \
		fi; \
	    fi; \
	done;
	$(buildAll)


ifeq (,$(findstring OSF,$(BFARCH)))
export INCPACKAGES := $(filter-out $(BUILDALLSKIP_include),$(PACKAGES))
export INCPACKAGES_1 := $(wordlist 1,300,$(INCPACKAGES))
export INCPACKAGES_2 := $(wordlist 301,600,$(INCPACKAGES))
export INCPACKAGES_3 := $(wordlist 601,900,$(INCPACKAGES))
export INCPACKAGES_4 := $(wordlist 901,$(words $(INCPACKAGES)),$(INCPACKAGES))
endif
export INCOPTIONS := NOPKG=yes
export MAKE

fastinclude:
	@echo "-> fastinclude:"
	packages=""; \
	for package in xx $(SCHEMALIST); do \
	    packages="$$package $$packages"; \
	done; \
	for package in $$packages; do \
	    if [ "$$package" != "xx" ]; then \
		$(MAKE) $$package.include NOPKG=yes; \
		if [ $$? -ne 0 -a "$(PRODUCTION_RELEASE)" != "" -a "$(SCHEMA_ERROR_IGNORE)" = "" ]; then \
		    exit 1; \
		fi; \
	    fi; \
	done;
	-@$(slaclog1) 
	makeParallelInc;
	-@$(slaclog2) 


#==========================================================================
#++ install directories.

# check real directory pointed to by a link created by 'newrel -s'
checklink:
	@if [ "$(DEBUG12)" ]; then set -x; fi ; \
	 for v in bin lib shlib tmp shtmp test results bin/$(BFARCH) lib/$(BFARCH) shlib/$(BFARCH) tmp/$(BFARCH) shtmp/$(BFARCH); do \
	    if [ -h $$v ]; then \
		if (ls -ltd $$v | grep '>' >/dev/null); then \
		    linkto=`ls -ldt $$v | awk '{print $$NF}' -`; \
		    if [ ! -d $$linkto ] ; then \
			echo "$$linkto no longer exist, remake it"; \
			mkdir -p $$linkto; \
		    fi; \
		fi; \
	    fi; \
	done

#++ winstalldirs: simply echo "-> installdirs:"
winstalldirs:
ifdef ONEDIR
	@[ -n "$(DEBUG23)" ] && echo "MAKEGOAL=/$(MAKEGOAL)/ argnum=$$#"; \
	[ "$(MAKEGOAL)" != "installdirs" ] && echo "-> $(PACKAGES).installdirs"; \
	echo "   Type 'gmake installdirs' if there is missing directory.";
else
	@[ -n "$(DEBUG23)" ] && echo "MAKEGOAL=/$(MAKEGOAL)/ argnum=$$#"; \
	echo "-> installdirs:"
endif

# create Sun's template as soft link to lib/<arch>/templates/<pkg>
ifneq (,$(findstring SunOS5,$(BFARCH)))
define suntemplates
		$(foreach var,$(PACKAGES),$(workdir)$(var)/SunWS_cache) \
		$(foreach var,$(PACKAGES),$(workdir)$(var)/obj/SunWS_cache) \
		$(foreach var,$(PACKAGES),$(workdir)$(var)/obj2/SunWS_cache) \
		$(foreach var,$(PACKAGES),$(shworkdir)$(var)/SunWS_cache) \
		$(foreach var,$(PACKAGES),$(shworkdir)$(var)/shobj/SunWS_cache) \
		$(foreach var,$(PACKAGES),$(shworkdir)$(var)/shobj2/SunWS_cache)
endef
endif

ifeq (,$(findstring NOInstalldirs,$(OPTIONS)))
installdirs: winstalldirs checklink $(bindir) $(shlibdir) $(libdir) $(libdir)templates $(shlibdir)templates $(shlibdir)ldlink $(shlibdir)LDLINK $(workdir) $(shworkdir) $(includedir) $(testdir) $(dbdir)/GNUmakefile results \
		$(foreach var,$(PACKAGES),$(workdir)$(var)) \
		$(foreach var,$(PACKAGES),$(workdir)$(var)/obj) \
		$(foreach var,$(PACKAGES),$(workdir)$(var)/obj2) \
		$(foreach var,$(PACKAGES),$(shworkdir)$(var)) \
		$(foreach var,$(PACKAGES),$(shworkdir)$(var)/shobj) \
		$(foreach var,$(PACKAGES),$(shworkdir)$(var)/shobj2) \
		$(suntemplates) \
		$(foreach var,$(PACKAGES),$(libdir)templates/$(var)) \
		$(foreach var,$(PACKAGES),$(shlibdir)templates/$(var)) \
		$(foreach var,$(PACKAGES),$(testdir)$(var)) \
		$(foreach var,$(PACKAGES),$(includedir)$(var))
else
installdirs:  
endif

$(dbdir)/GNUmakefile: $(BFDIST)/releases/$(BFCURRENT)/SoftRelTools/GNUmakefile.database
	if [ ! -d $(@D) ] ; then \
	    mkdir $(@D); \
	fi
	if [ -r $< ]; then \
	    if [ -r $@ ]; then \
		rm -f $@; \
	    fi; \
	    echo "Creating database/GNUmakefile from release $(BFCURRENT)"; \
	    cp -f $< $@; \
	else \
	    if [ -r $@ ]; then \
		echo "Unable to update database/GNUmakefile from release $(BFCURRENT)"; \
	    else \
		echo "Unable to create database/GNUmakefile from release $(BFCURRENT)"; \
	    fi; \
	    if [ -d SoftRelTools ]; then \
		echo "Updating database/GNUmakefile from local SoftRelTools package"; \
		cp -f SoftRelTools/GNUmakefile.database $@; \
	    fi; \
	fi

$(BFDIST)/releases/$(BFCURRENT)/SoftRelTools/GNUmakefile.database:
	echo "$@ doesn't exist!"

results:
	if [ ! -d results ] ; then mkdir results ; fi

checkdep:
	@echo "-> checkdep: looking for packages to recompile:"
	checkdep

.PHONY:	installdirs

#__________________________________________________________________________
#++ target to define directories

$(bindir): 
	if [ ! -d bin ] ; then mkdir bin ; fi
	- mkdir $(bindir)

$(shlibdir):
	if [ ! -d shlib ] ; then mkdir shlib ; fi
	- mkdir $(shlibdir)

$(shlibdir)tmp:
	if [ ! -d shlib/tmp ] ; then mkdir shlib/tmp ; fi
	- mkdir $(shlibdir)tmp

$(libdir): 
	if [ ! -d lib ] ; then mkdir lib ; fi
	- mkdir $(libdir)

$(libdir)templates:
	- mkdir $(libdir)templates

$(shlibdir)templates:
	- mkdir $(shlibdir)templates

$(shlibdir)ldlink:
	- mkdir $(shlibdir)ldlink

$(shlibdir)LDLINK:
	- mkdir $(shlibdir)LDLINK

define mktmp
	if [ "$(DEBUG34)" ]; then set -x; fi; \
	remark="remove bad link"; rm $$wkdir 2>/dev/null; \
	topdir=$(notdir $(TOPDIR)); \
	scratch=`cat .scratch 2>/dev/null`; \
	if [ -z "$$scratch" ]; then \
	    mkdir $$wkdir; \
	else \
	    linkto=1; \
	    if [ "$(BFSITE)" = "slac" ]; then \
	        ans="`cd $$scratch; /bin/pwd`"; \
		if [ `expr "$$ans" : "/a/.*" ` -gt 0 ]; then \
		   if [ `expr "$(BFARCH)" : '.*SunOS5.*'` -gt 0 ]; then \
		       if [ "$$work" = "tmp" ]; then echo "leave tmp and shtmp on the test release for SunOS"; fi; \
		       if [ ! -d $$wkdir ]; then mkdir $$wkdir; fi; \
		       linkto=; \
		   fi; \
                fi; \
	    fi; \
	    if [ "$$linkto" ]; then \
	        if [ ! -d $$scratch/$$topdir/$$work/$(BFARCH) ]; then mkdir -p $$scratch/$$topdir/$$work/$(BFARCH); fi; \
		if [ ! -d $$work ]; then mkdir $$work; fi; \
		ln -s $$scratch/$$topdir/$$work/$(BFARCH) $$work/$(BFARCH); \
	    fi; \
	fi
endef

$(workdir): 
	if [ ! -d tmp ] ; then mkdir tmp ; fi
	-work=tmp; wkdir=$(workdir); $(mktmp)

$(shworkdir): 
	if [ ! -d shtmp ] ; then mkdir shtmp ; fi
	-work=shtmp; wkdir=$(shworkdir); $(mktmp)
#	- mkdir $(shworkdir)

$(includedir): 
	if [ ! -d include ] ; then mkdir include ; fi

$(foreach var,$(PACKAGES),$(workdir)$(var)):
	- mkdir $@

$(foreach var,$(PACKAGES),$(workdir)$(var)/obj):
	- mkdir $@

$(foreach var,$(PACKAGES),$(workdir)$(var)/obj2):
	- mkdir $@

$(foreach var,$(PACKAGES),$(shworkdir)$(var)/shobj):
	- mkdir $@

$(foreach var,$(PACKAGES),$(shworkdir)$(var)/shobj2):
	- mkdir $@

ifneq (,$(findstring SunOS5,$(BFARCH)))
$(foreach var,$(PACKAGES),$(workdir)$(var)/SunWS_cache):
	- if [ ! -r "$@" ]; then ln -s $(libdir)templates/$(notdir $(@D)) $@; fi

$(foreach var,$(PACKAGES),$(workdir)$(var)/obj/SunWS_cache):
	- if [ ! -r "$@" ]; then ln -s $(libdir)templates/$(word $(words $(subst /, ,$(dir $(@D)))), $(subst /, ,$(dir $(@D)))) $@; fi

$(foreach var,$(PACKAGES),$(workdir)$(var)/obj2/SunWS_cache):
	- if [ ! -r "$@" ]; then ln -s $(libdir)templates/$(word $(words $(subst /, ,$(dir $(@D)))), $(subst /, ,$(dir $(@D)))) $@; fi

$(foreach var,$(PACKAGES),$(shworkdir)$(var)/SunWS_cache):
	- if [ ! -r "$@" ]; then ln -s $(shlibdir)templates/$(notdir $(@D)) $@; fi

$(foreach var,$(PACKAGES),$(shworkdir)$(var)/shobj/SunWS_cache):
	- if [ ! -r "$@" ]; then ln -s $(shlibdir)templates/$(word $(words $(subst /, ,$(dir $(@D)))), $(subst /, ,$(dir $(@D)))) $@; fi

$(foreach var,$(PACKAGES),$(shworkdir)$(var)/shobj2/SunWS_cache):
	- if [ ! -r "$@" ]; then ln -s $(shlibdir)templates/$(word $(words $(subst /, ,$(dir $(@D)))), $(subst /, ,$(dir $(@D)))) $@; fi
endif

$(foreach var,$(PACKAGES),$(libdir)templates/$(var)):
	- mkdir $@

$(foreach var,$(PACKAGES),$(shlibdir)templates/$(var)):
	- mkdir $@

$(foreach var,$(PACKAGES),$(shworkdir)$(var)):
	- mkdir $@

$(testdir):
	if [ ! -d test ] ; then mkdir test ; fi
	- mkdir $(testdir)

$(foreach var,$(PACKAGES),$(testdir)$(var)):
	- mkdir $@

$(foreach var,$(PACKAGES),$(includedir)$(var)):
	-pkg=$(notdir $@); \
	if [ "$$pkg" = "database" ]; then \
	    remark="skip database"; \
	elif [ -d $$pkg/$$pkg ]; then \
	    ln -s ../$$pkg/$$pkg $@; \
	elif [ -d $$pkg/include ]; then \
	    ln -s ../$$pkg/include $@; \
	elif [ -d $$pkg/src ]; then \
	    ln -s ../$$pkg/src $@; \
	else \
	    ln -s ../$$pkg/ $@; \
	fi;

.SUFFIXES = 

#==========================================================================
#++ flags to be passed over to the sub-make.

ifdef BFSITE
	siteconf=-I$(BFDIST)/releases/$(BFCURRENT)/SiteConfig/$(BFSITE)
	localconf= -I$(TOPDIR)/SiteConfig/$(BFSITE)
endif

#-> search path of Makefiles, etc
export OVERRIDES = shlibdir=$(shlibdir) libdir=$(libdir) bindir=$(bindir) LDFLAGS="$(LDFLAGS)" workdir=$(workdir)$(basename $@)/ shworkdir=$(shworkdir)$(basename $@)/ templatedir=$(templatedir)$(basename $@)/ testdir=$(testdir)$(basename $@)/ -I$(BFOVERRIDE)/$(BFCURRENT) $(localconf) $(siteconf) -I$(TOPDIR)/include -I$(BFDIST)/releases/$(BFCURRENT)/include mandir=$(mandir) docdir=$(docdir) CPPFLAGS="$(CPPFLAGS)" IDLINCLFLAGS="$(IDLINCLFLAGS)"

# this null item is needed to avoid passing it down each time
GNUmakefile:

#-> function to display current time 
ifneq (,$(findstring Time,$(OPTIONS)))
define qtime
(`date; uptime`)
endef
else
define qtime
(`date`)
endef
endif

#-> function to display elapsed time 
ifneq (,$(findstring Time,$(OPTIONS)))
define etime
$(ETIME)
endef

#-> function to display end of target
define endtime
echo "-> end of $@"
endef
endif

#==========================================================================
#++ x.y style targets.

# Make sure that the database pseudo-package & GNUmakefile are 
# automatically created if necessary.
database.%: $(dbdir)/GNUmakefile
	@echo "-> $@:   $(qtime)"; pkg=$(word 1,$(subst ., ,$@)); $(MAKE) -C $(notdir $(basename $@)) $(OVERRIDES) $(word 2,$(subst ., ,$@)) MAKEGOAL=$(word 2,$(subst ., ,$@)) srcdir=$(shell /bin/pwd)/$$pkg PACKAGE=$$pkg

database@%: $(dbdir)/GNUmakefile
	@echo "-> $@:   $(qtime)"; pkg=$(word 1,$(subst ., ,$@)); $(MAKE) -C $(notdir $(basename $(subst @,.,$@))) $(OVERRIDES) $(word 2,$(subst @, ,$@)) MAKEGOAL=$(word 2,$(subst @, ,$@)) srcdir=$(shell /bin/pwd)/$$pkg PACKAGE=$$pkg

# anything of the form <package>@<item> means to do a gmake item
# in the package directory - this is a newer form, allowing periods in
# <item> and more than one <package> in the target
$(foreach var,$(PACKAGES),$(var)@%) : 
	@pkg=$(word 1,$(subst ., ,$@)); $(MAKE) -C $(notdir $(basename $(subst @,.,$@))) $(OVERRIDES) $(word 2,$(subst @, ,$@)) MAKEGOAL=$(word 2,$(subst @, ,$@)) srcdir=$(shell /bin/pwd)/$$pkg PACKAGE=$$pkg

#__________________________________________________________________________
#++ the <pkg>.<target> style targets. 
#   anything of the form <package>.<item> means to do a gmake item
#   in the package directory.

$(foreach var,$(PACKAGES),$(var).%) : 
	-@$(slaclog1)
	@echo "-> $@:   $(qtime)"; pkg=$(word 1,$(subst ., ,$@)); goal=$(word 2,$(subst ., ,$@)); if [ "$(SUBDIR)" ]; then Subdir="/$(SUBDIR)"; fi; $(etime) $(MAKE) -C $(notdir $(basename $@))$$Subdir $(OVERRIDES) $$goal MAKEGOAL=$$goal srcdir=$(shell /bin/pwd)/$$pkg PACKAGE=$$pkg; $(endtime)
	-@$(slaclog2)

# any target not recognized gets passed to all subdirectories
# i.e.  gmake foo  means do "gmake foo" in all subdirectories
% :
	-@$(slaclog1)
	@pkg=$(word 1,$(subst ., ,$@)); \
	if [ "$$pkg" = "" ]; then \
	    for package in xx $(PACKAGES); do \
		if [ "$$package" != "xx" ]; then \
		    echo "-> $$package.$@"; \
		    $(MAKE) srcdir=$(shell /bin/pwd) $$package.$@ NOPKG=yes; \
		fi; \
	    done; \
	fi
	-@$(slaclog2)

#__________________________________________________________________________
#-> SRT spec function 
spec specver specpath:
	@echo "-> $@"; \
	if [ "$$specmk" = "" ]; then \
		echo "can not find spec.mk, quit" >&2; \
		exit 2; \
	fi; \
	if [ "$(DEBUG)" != "" ]; then \
		echo "use $$specmk"; \
	fi; \
	$(MAKE) -C /tmp -f $$specmk $(OVERRIDES) $@

#==========================================================================
#++ one time only section. Set up SRT environment.
#   This section can be positioned near the end of this file.

ifeq (,$(NOTFIRST))
ifdef MAKEOVERRIDES
    override MAKEOVERRIDES += NOTFIRST=1
else
    export NOTFIRST := 1
endif

    ERROR := $(shell $(MAKE) --version | head -1 | cut -f1-4 -d" " >& 2; \
	       echo "Build OPTIONS = $(OPTIONS)" >& 2)

ifdef DEBUG02
    ERROR := $(shell echo SCHEMAPACK=$(SCHEMAPACK) >& 2)
    ERROR := $(shell echo SCHEMALIST=$(SCHEMALIST) >& 2)
endif   ## DEBUG02

#__________________________________________________________________________
#++ temporary patch, warn user if they have -native/-Objy in BFARCH

ifneq (,$(findstring native,$(BFARCH)))
   ERROR := $(shell echo "Warning: you have -native in BFARCH, hope it's all right" >&2)
endif
ifneq (,$(findstring Objy,$(BFARCH)))
   ERROR := $(shell echo "Warning: you have -Objy in BFARCH, hope it's all right" >&2)
endif

#__________________________________________________________________________
#++ get architecture with 'uname -a'

    ERROR := $(shell echo "`uname -a`  [uname -a]" >& 2)

#__________________________________________________________________________
#++ Check bin/$BFARCH for installdirs

ifeq (,$(wildcard ./bin/$(BFARCH)))
ifeq (,$(findstring installdirs,$(MAKECMDGOALS)))
    ERROR := $(shell echo "Warning: bin/$(BFARCH) not found, may need to run gmake installdirs  [check-1]" >& 2)
endif	# not installdirs  
endif   # no bin/$(BFARCH)

#__________________________________________________________________________
#++ check PATH

     ifeq (,$(findstring ./bin/$(BFARCH),$(PATH)))
        srtpatherr:=1	
     endif
     ifeq (,$(findstring $(BFDIST)/releases/$(BFCURRENT)/bin/$(BFARCH),$(PATH)))
        srtpatherr:=1	
     endif
    ifdef srtpatherr
        ERROR:= $(shell echo "[Warning]: ./bin/$(BFARCH) (and/or) $(BFDIST)/releases/$(BFCURRENT)/bin/$(BFARCH) is not in PATH, type 'srtpath' to fix PATH." >& 2)
    endif	 

#__________________________________________________________________________
#++ check platform for incremental linking

    ifneq (,$(INCR))
        ifeq (,$(findstring AIX,$(BFARCH)))
            ERROR := $(shell echo Note: Incremental linking is supported only on AIX for now. INCR ignored. >& 2) 
            override export INCR :=
        endif
    endif

#__________________________________________________________________________
#++ check basic environment variables

    ifndef BFDIST
       ERROR := $(shell echo 'Error: BFDIST undefined' >& 2)
       exit
    endif

    ifndef BFARCH
       ERROR := $(shell echo 'Error: BFARCH undefined' >& 2)
       exit
    endif

    ifndef BFSITE
       ERROR := $(shell echo 'Note: BFSITE undefined, no site configuration will be used' >& 2)
    endif

#__________________________________________________________________________
#++ Check ldlink environment

ifneq (,$(findstring -Ldlink,$(OPTIONS)))
    ifeq (,$(wildcard $(BFDIST)/releases/$(BFCURRENT)/shlib/$(BFARCH)/ldlink/lib*_$(BFCURRENT).$(SOEXT)))
       item:=(1) 
       ldlinkmsg:= (2) $(BFDIST)/releases/$(BFCURRENT)/shlib/$(BFARCH)/ldlink is empty.
    endif

    ifeq (,$(findstring ./shlib/$(BFARCH)/ldlink,$(LD_LIBRARY_PATH)))
        ldlinkerr:=1	
    endif
    ifeq (,$(findstring $(BFDIST)/releases/$(BFCURRENT)/shlib/$(BFARCH)/ldlink,$(LD_LIBRARY_PATH)))
        ldlinkerr:=1
    endif
    ifneq (, $(ldlinkerr)$(ldlinkmsg))
#        ERROR:= $(shell echo "Warning: $(item) ./bin/$(BFARCH)/ldlink (and/or) $(BFDIST)/releases/$(BFCURRENT)/shlib/$(BFARCH)/ldlink is not in LD_LIBRARY_PATH, type 'srtpath' to fix LD_LIBRARY_PATH. $(ldlinkmsg)" >& 2)
    endif	 
endif

#__________________________________________________________________________
#++ give warning if third party software is defined in personal environment.
#   using site configuration is preferred.

ifneq (,$(OBJYBASE))
    ERROR := $(shell echo Note: "OBJYBASE is pre-defined in your personal environment as $(OBJYBASE)." >&2)
    ERROR := $(shell echo "      Make sure it is correct for this site." >&2)
endif

#__________________________________________________________________________
#++ include all mk files to be used on the top-level.
#   macros defined this way are global.

PKG = SoftRelTools
MKFILES = profile.mk profile.mk+
-include $(foreach MK, $(MKFILES), $(shell $(findmk)))
MK=spec.mk
export specmk:=$(shell $(findmk))
endif
