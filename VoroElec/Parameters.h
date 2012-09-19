/*
 * Parameters.h
 *
 *  Created on: Mar 12, 2012
 *      Author: marchi
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#include <string>
#include <cmath>

namespace Parameters{
	struct Running {
		static int nframes;
	};
	struct Input{
		static float dx;
		static float cut;
		static int Dir;
		static std::string filename;
		static bool verbose;
		static bool PrintRoundIon;
		static float freq;
		static bool bRot;
		static bool bTiming;
		static bool bVoronoi;
		static bool bTest;
		static bool bTestConv;
		static int nAverage;
	};
	const double pi=3.14159265358979323846,
			twopi=2.0*pi,
			avogad=6.0225e23,boltz=1.38054e-23,
			gascon=8.3143,
    		planck=6.6256e-34,
    		elechg=1.602176487e-19,
    		epso=8.85418782e-12,
    		kel=1.0/(epso*4.0*pi),
    		boxl=2.0,
    		unitm=1.0/(avogad*1000.0),
    		unitl=1.0e-9,
    		unitt=1.0e-15,
    		unite=unitm*(unitl/unitt)*(unitl/unitt),
    		unitc=4.0*pi*epso*unitl*unite/(elechg*elechg),
    		unitp=(unite/unitl*unitl*unitl)/1.0e6,
    		efact=unite*avogad,
    		efact_nm=100*unite*avogad,
    		hartree=4.35981 * 1.0e-18,
    		unitepot=unite/sqrt(unitc)/hartree,
    		lbohr=0.52917706,
    		unitefield=unitepot*lbohr,
    		unitfield=kel*elechg/unitl,
    		kT300=gascon*300.0/1000.0;
}




#endif /* PARAMETERS_H_ */
