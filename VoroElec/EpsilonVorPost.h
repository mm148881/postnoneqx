/*
 * EpsilonVorPost.h
 *
 *  Created on: Mar 30, 2012
 *      Author: marchi
 */

#ifndef EPSILONVORPOST_H_
#define EPSILONVORPOST_H_
#include "EpsilonVor.h"
#include "Parameters.h"
#include <iostream>
#include <fstream>
//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/uniform_int_distribution.hpp>

namespace EpsilonNS {
//static boost::mt19937 gen;
static int roll_dice(int i, int j) {
//    boost::random::uniform_int_distribution<> dist(i, j);
//    return dist(gen);
	return 0;
}

class EpsilonVorPost: public EpsilonVor {
protected:
	vector<IndexedMatrix> Data;
	virtual void WriteIt(std::ofstream &) const;
public:
	EpsilonVorPost();
	virtual void Rdf(double,double);
	virtual string getLabel(int i){return ResLabel[i];}
	virtual void Statistics();
	virtual ~EpsilonVorPost();
};

} /* namespace EpsilonNS */
#endif /* EPSILONVORPOST_H_ */
