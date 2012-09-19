/*
 * EpsilonGridPost.h
 *
 *  Created on: Mar 30, 2012
 *      Author: marchi
 */

#ifndef EPSILONGRIDPOST_H_
#define EPSILONGRIDPOST_H_
#include "EpsilonGrid.h"
#include "Parameters.h"
#include <iostream>
#include <fstream>

namespace EpsilonNS {
class EpsilonGridPost: public EpsilonGrid {
protected:
	vector<IndexedMatrix> RdfData;
	virtual void WriteIt(std::ofstream &) const;
public:
	EpsilonGridPost();
	virtual void Rdf(const double , const double );
	virtual ~EpsilonGridPost();
};

} /* namespace EpsilonNS */
#endif /* EPSILONGRIDPOST_H_ */
