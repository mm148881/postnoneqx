# Makefile for 
# Copyright (C) 2012 Massimo Marchi <Massimo.Marchi@cea.fr>
# $Header$
# $Log$
# what flags you want to pass to the C compiler & linker
AM_CXXFLAGS = -g0
AM_LDFLAGS = -L$(top_builddir)/alglib -L$(top_builddir)/MyLibrary -L$(top_builddir)/VoroElec
AM_CPPFLAGS = -I$(top_builddir)/alglib -I$(top_builddir)/math.nist/jama -I$(top_builddir)/math.nist/tnt1 -I$(top_builddir)/MyLibrary -I$(top_builddir)/VoroElec
LDADD = -lmy -lVoroElec -lalglib 

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = PostNonEq
PostNonEq_SOURCES = EpsilonGridNE.cpp EpsilonGridNE.h EpsilonVorNE.cpp EpsilonVorNE.h\
	PostNonEq.cpp

