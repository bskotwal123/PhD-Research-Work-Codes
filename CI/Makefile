#-----------------------------------------------------------------------------
# Description: Makefile to build libCI.so
# Created:     Oct 2014
# Author:      Shakepeare's ghost
#-----------------------------------------------------------------------------
ifndef ROOTSYS
$(error *** Please set up Root)
endif

name    := CI

# Sub-directories
srcdir	:= src
libdir	:= lib
incdir	:= include

$(shell mkdir -p lib)

# Set this equal to the @ symbol to suppress display of instructions
# while make executes
ifdef verbose
AT 	:=
else
AT	:= @
endif

srcs	:= 	$(srcdir)/JetCorrectionUncertainty.cc \
		$(srcdir)/JetCorrectorParameters.cc \
		$(srcdir)/SimpleJetCorrectionUncertainty.cc \
		$(srcdir)/CIXsection.cc \
		$(srcdir)/CISpectrum.cc \
		$(srcdir)/QCDSpectrum.cc \
		$(srcdir)/JetSpectrum.cc \
		$(srcdir)/JetSpectrumSmeared.cc \
		$(srcdir)/hutil.cc \
		$(srcdir)/Bayes.cc \
		$(srcdir)/PDFWrapper.cc \
		$(srcdir)/PDFunction.cc \
		$(srcdir)/RooMultinomial.cc \
		$(srcdir)/RooInclusiveJetPdf.cc

othersrcs	:= $(srcdir)/Utilities.cc

dictsrcs:= $(srcs:.cc=_dict.cc)
sources	:= $(srcs) $(othersrcs) $(dictsrcs)
objects	:= $(sources:.cc=.o)

sharedlib := $(libdir)/lib$(name).so

# Display list of applications to be built
#say	:= $(shell echo -e "Apps: $(applications)" >& 2)
#say	:= $(shell echo -e "AppObjs: $(appobjs)" >& 2)
#say	:= $(shell echo -e "Objects: $(objects)" >& 2)
#$(error bye!) 

#-----------------------------------------------------------------------
# 	Define which compilers and linkers to use

# 	C++ Compiler
CXX	:= clang++

CINT    := rootcint

# 	Define paths to be searched for C++ header files (#include ....)
CPPFLAGS:= -I. -I$(incdir) -I$(srcdir) $(shell root-config --cflags)

# 	Define compiler flags to be used
#	-c		perform compilation step only 
#	-g		include debug information in the executable file
#	-O2		optimize
#	-ansi	require strict adherance to C++ standard
#	-Wall	warn if source uses any non-standard C++
#	-pipe	communicate via different stages of compilation
#			using pipes rather than temporary files

CXXFLAGS:= -c -g -O2 -ansi -Wall -pipe -fPIC

#	C++ Linker
#   set default path to shared library
LD	:= clang++

OS	:= $(shell uname -s)
ifeq ($(OS),Darwin)
	LDSHARED	:= $(LD) -dynamiclib
else
	LDSHARED	:= $(LD) -shared
endif

#	Linker flags
LDFLAGS := -g $(root-config --ldflags)

# 	Libraries

LIBS	:=  $(shell root-config --libs) \
-L$(libdir) -lRooStats -lRooFit -lRooFitCore -lPyROOT \
-lMinuit  -lMathMore -lMathCore


#	Rules
#	The structure of a rule is
#	target : source
#		command
#	The command makes a target from the source. 
#	$@ refers to the target
#	$< refers to the source

lib:	$(sharedlib)

# Syntax:
# list of targets : target pattern : source pattern

$(sharedlib)	: $(objects)
	@echo ""
	@echo "---> Linking `basename $@`"
	$(AT)$(LDSHARED) $(LDFLAGS) -fPIC $(objects) $(LIBS) -o $@

$(objects)	: %.o	: %.cc
	@echo ""
	@echo "---> Compiling `basename $<`" 
	$(AT)$(CXX) $(CXXFLAGS) $(CPPFLAGS)  $< -o $@ >& $*.FAILED
	$(AT)rm -rf $*.FAILED

$(dictsrcs)	: $(srcdir)/%_dict.cc : $(incdir)/%.h
	@echo ""
	@echo "---> Building dictionary from `basename $^`" 
	$(AT)$(CINT) -f $@ -c $^

# 	Define clean up rules
clean   :
	rm -rf $(libdir)/* $(srcdir)/*_dict*.* $(srcdir)/*.o
