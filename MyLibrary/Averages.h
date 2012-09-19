/*
 * Averages.h
 *
 *  Created on: Jun 21, 2011
 *      Author: marchi
 */

#ifndef AVERAGES_H_
#define AVERAGES_H_
#include "RhoAverage.h"
#include "RhoDipAverage.h"

class Averages {
public:
	RhoAverage ** Ro;
	RhoDipAverage * RoDip;
	RhoAverage * RhoSqrK;
	RhoAverage * RhoSqrK_Dip;
	RhoAverage * RhoSqrK_Ratio;
	Averages();
	virtual ~Averages();
};

#endif /* AVERAGES_H_ */
