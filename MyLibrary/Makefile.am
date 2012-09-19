# Makefile for 
# Copyright (C) 2012 Massimo Marchi <Massimo.Marchi@cea.fr>
# $Header$
# $Log$
# what flags you want to pass to the C compiler & linker
AM_CXXFLAGS = -g0
AM_LDFLAGS = 
AM_CPPFLAGS = -I$(top_builddir)/alglib -I$(top_builddir)/math.nist/jama -I$(top_builddir)/math.nist/tnt1 -I$(top_builddir)/VoroElec

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
lib_LIBRARIES = libmy.a
libmy_a_SOURCES = Array.h AtomIndex.cpp AtomIndex.h Atoms.cpp Atoms.h\
	AtomsPDB.cpp AtomsPDB.h Averages.cpp Averages.h BSplineCoeffs.cpp\
	BSplineCoeffs.h BSpline.cpp BSpline.h BspModules.cpp BspModules.h\
	CellSymmetry.cpp CellSymmetry.h Charges.cpp Charges.h Density.cpp\
	Density.h Dielectric.h DiffCoeffs.cpp DiffCoeffs.h DiffCoeffs_Spec.hpp\
	Dipoles.cpp Dipoles.h EpsofK.cpp EpsofK.h fftw++.cpp\
	fftw++.h ForceConst.cpp ForceConst.h gmx_polar.cpp gmx_polar.h\
	Grid.cpp Grid.h Grid_SpecTemplates.hpp Metric.cpp Metric.h\
	MyFiles.h Mymap.cpp Mymap.h Mymap_Spec.hpp MyStrings.cpp\
	MyStrings.h MyVec.h NeighList.cpp NeighList.h\
	Phi.cpp Phi.h ResidueAvg.cpp ResidueAvg.h ResidueCM.cpp\
	ResidueCM.h RhoAverage.cpp RhoAverage.h RhoAvg.cpp RhoAvg.h\
	Rho.cpp RhoDipAverage.cpp RhoDipAverage.h RhoDip.cpp RhoDip.h\
	Rho.h simpson.hpp WeightedDensity.cpp WeightedDensity.h
