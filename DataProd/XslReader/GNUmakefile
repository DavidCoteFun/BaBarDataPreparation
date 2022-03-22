# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile,v 1.2 2004/10/05 23:41:39 cote Exp $
#################################################################
#++ library products				[build it with 'lib']

LIBREMOVEFILES := 
LIBTMPLFILES := 
LIBDDLORDERED := 

#################################################################
#++ extra binary products	[not in production, build it with extrabin]

EXTRABINS := BetaMicroApp
$(addprefix $(bindir),$(EXTRABINS)): $(bindir)% : AppMicroBuild.o

#################################################################
#++ binary products				[build it with 'bin']

BINS := BetaMiniApp

BINCCFILES := AppUserBuild.cc

#++ Binary rules		 [in production, build it with 'bin']

$(bindir)BetaMiniApp: AppUserBuild.o
$(bindir)BetaMicroApp: AppMicroBuild.o

#++ shell script products.. 			[build it with 'bin']
BINSCRIPTS := 

#################################################################
#++ regression test scripts			[build it with 'test']


#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk

