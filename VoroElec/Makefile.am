# Makefile for 
# Copyright (C) 2012 Massimo Marchi <Massimo.Marchi@cea.fr>
# $Header$
# $Log$
# what flags you want to pass to the C compiler & linker
AM_CXXFLAGS = -g0
AM_LDFLAGS = -L$(top_builddir)/alglib 
AM_CPPFLAGS = -I$(top_builddir)/alglib -I$(top_builddir)/math.nist/jama -I$(top_builddir)/math.nist/tnt1 -I$(top_builddir)/MyLibrary
LIBS  = -lalglib

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
lib_LIBRARIES = libVoroElec.a
libVoroElec_a_SOURCES = Accumulate.cpp Accumulate.h Analyze_b.h ChangeInput.cpp\
	ChangeInput.h Communicator.hpp Epsilon.cpp EpsilonGrid.cpp EpsilonGrid.h\
	EpsilonGridPost.cpp EpsilonGridPost.h Epsilon.h EpsilonVor.cpp EpsilonVor.h\
	EpsilonVorPost.cpp EpsilonVorPost.h Euler.cpp Euler.h Field.cpp\
	Field.h Flag.cpp Flag.h gmx_voroelec.h gmx_voronoi.h\
	GridEps.cpp GridEps.h Gridn.cpp Gridn.h Gridn_SpecTemplates.hpp\
	LocalEps.cpp LocalEps.h MyDipole.cpp MyDipole.h MyMPI.cpp\
	MyMPI.h MyMPI_Spec.hpp MyUtilClass.cpp MyUtilClass.h Parameters.cpp\
	Parameters.h Polarization.cpp Polarization.h Residue.cpp Residue.h\
	Voronoi.cpp Voronoi.h VoronoiPrint.cpp VoronoiPrint.h
