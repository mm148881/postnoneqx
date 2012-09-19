/*
 * EpsilonGridNE.h
 *
 *  Created on: Jun 28, 2012
 *      Author: marchi
 */

#ifndef EPSILONGRIDNE_H_
#define EPSILONGRIDNE_H_

#include "EpsilonGridPost.h"

namespace EpsilonNS {
class EpsilonGridNE: public EpsilonGridPost {
	int Dir;
	virtual void WriteIt(std::ofstream &) const;
public:
	EpsilonGridNE();
	virtual void Rdf(const double , const double );
	virtual void setdir(const int D){Dir=D;};
	virtual ~EpsilonGridNE();
};

} /* namespace EpsilonNS */
#endif /* EPSILONGRIDNE_H_ */
