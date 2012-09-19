/*
 * RhoDip.h
 *
 *  Created on: Jun 21, 2011
 *      Author: marchi
 */

#ifndef RHODIP_H_
#define RHODIP_H_

#include "Grid.h"
#include "Dipoles.h"
#include "typedefs.h"
#include <typeinfo>
#include <cstring>

class RhoDip: public Grid<DIM> {
public:
	RhoDip():Grid<DIM>::Grid(){};
	RhoDip & operator=(const Grid<DIM> &);
	RhoDip & operator=(const double &);
	void Density(Dipoles & );
	virtual ~RhoDip(){};
	virtual void Accumulate(){};

	friend std::ofstream & operator<<(std::ofstream &, RhoDip &);
	friend std::ifstream & operator>>(std::ifstream &, RhoDip & );
};

#endif /* RHODIP_H_ */
