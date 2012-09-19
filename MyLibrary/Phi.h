/*
 * Phi.h
 *
 *  Created on: May 24, 2011
 *      Author: marchi
 */

#ifndef PHI_H_
#define PHI_H_
#include <iostream>
#include <typedefs.h>
#include "Rho.h"
#include "fftw++.h"
#include "Grid.h"
using namespace Array;
using namespace fftwpp;
using std::cout;
using std::endl;


class Phi : public Grid<1> {
	void GetPhi(const Rho &);
public:
	Phi(){};
	Phi(const Rho &);
	Phi & operator=(const Grid<1> &);
	Phi & operator=(array3<Complex> &);
	Phi & operator=(array3<double> &);
	Phi & operator=(const double & );
	Phi & operator()(const Rho &);
};

#endif /* PHI_H_ */
