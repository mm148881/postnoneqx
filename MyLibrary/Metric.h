/*
 * Metric.h
 *
 *  Created on: May 25, 2011
 *      Author: marchi
 */

#ifndef METRIC_H_
#define METRIC_H_
#include "typedefs.h"

class Metric {
	matrix co, oc;
	static const matrix Idtity;
	void invertCO();
	int count;
public:
	Metric();
	Metric(const matrix &);
	Metric(const Metric &);
	virtual ~Metric();
	Metric & operator()(const matrix &);
	Metric & operator()(const Metric &);
	Metric & operator=(const Metric &);
	Metric & operator+=(const Metric &);
	Metric operator/(const double);
	const matrix & getCO() const {return co;};
	const matrix & getOC() const {return oc;};

	double getVol(){ return (double) (co[XX][XX]*(co[YY][YY]*co[ZZ][ZZ]-co[ZZ][YY]*co[YY][ZZ])
		  -co[YY][XX]*(co[XX][YY]*co[ZZ][ZZ]-co[ZZ][YY]*co[XX][ZZ])
		  +co[ZZ][XX]*(co[XX][YY]*co[YY][ZZ]-co[YY][YY]*co[XX][ZZ]));};
};

#endif /* METRIC_H_ */
