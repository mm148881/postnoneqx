/*
 * ResidueAvg.h
 *
 *  Created on: Sep 21, 2011
 *      Author: marchi
 */

#ifndef RESIDUEAVG_H_
#define RESIDUEAVG_H_
#include "typedefs.h"
#include "ResidueCM.h"
#include "Array.h"
#include <iostream>
using namespace Array;
using std::cout;
using std::endl;

class ResidueAvg {
protected:
	unsigned int nr;
	int nr0;
	Array2<double> Avg;
	Array2<double> Coord;
	Array4<double> Avg2;
	Array4<double> Rms;
	static int counter;
public:
	ResidueAvg();
	ResidueAvg(int);
	ResidueAvg & operator+=(ResidueCM &);
	void Spit();
	void Average();
	virtual ~ResidueAvg();
	Array4<double> & getRms(){return Rms;};
	Array2<double> & getCoord(){return Coord;};
};

#endif /* RESIDUEAVG_H_ */
