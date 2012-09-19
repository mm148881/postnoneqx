/*
 * gmx_voroelec.h
 *
 *  Created on: Dec 27, 2011
 *      Author: marchi
 */

#ifndef GMX_VOROELEC_H_
#define GMX_VOROELEC_H_

#include "MyDipole.h"
#include "Polarization.h"
#include "GridEps.h"
#include <string>
#include "RhoAverage.h"
#include "EpsilonGrid.h"
#include "EpsilonVor.h"

using namespace PolarizationNS;
using namespace GridEpsNS;

namespace gmx_voroelec {
	std::ofstream fout;
	MyDipole * Dips;
	Polarization * P;
	Polarization * P_t;
	Polarization * Pa_t;
	Field * E;
	Field * E_t;
	Field * Ea_t;
	Rho * Rho_t;
	RhoAverage * Rhoa_t;
	EpsilonNS::Epsilon * MyEps0;
	static ofstream opfout;
	static vector<double> * dataIn;
	static vector<double> * dataOut;
}

#endif /* GMX_VOROELEC_H_ */
