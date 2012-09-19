/*
 * Density.cpp
 *
 *  Created on: Oct 3, 2011
 *      Author: marchi
 */

#include "Density.h"
Density & Density::operator=(const Grid<1> & x){
	Grid<1> & My=*this;
	My=x;
	return *this;
}
Density & Density::operator=(const double & x){
	Grid<1> & My=*this;
	My=x;
	return *this;
}

Density::~Density() {
	// TODO Auto-generated destructor stub
}

