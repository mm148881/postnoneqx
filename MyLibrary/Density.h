/*
 * Density.h
 *
 *  Created on: Oct 3, 2011
 *      Author: marchi
 */

#ifndef DENSITY_H_
#define DENSITY_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "typedefs.h"
#include "Grid.h"
#include <typeinfo>
#include <cstring>
using namespace Array;

class Density : public Grid<1> {
protected:
public:
	Density():Grid<1>::Grid(){};
	Density & operator=(const Grid<1> &);
	Density & operator=(const double &);
	double Sum(){
		double tot=0.0;
		for(unsigned int i=0;i<this->Size();i++) tot+=this->v[i];
		return tot;
	};
	virtual ~Density();
};

#endif /* DENSITY_H_ */
