/*
 * ForceConst.h
 *
 *  Created on: Sep 22, 2011
 *      Author: marchi
 */

#ifndef FORCECONST_H_
#define FORCECONST_H_
#include "jama_svd.h"
#include "Array.h"
#include "ResidueAvg.h"
#include "typedefs.h"
#include <iostream>
#include <fstream>


using namespace JAMA;
using namespace TNT;
using namespace Array;

class ForceConst {
protected:
	Array2D<double> Rms;
	Array2D<double> Dist;
	Array2D<double> F;
	Array2D<double> transpose(Array2D<double>);
	float Percent;
public:
	ForceConst();
	ForceConst(ResidueAvg &);
	double * operator[](int n){return F[n];};
	void Svd();
	void SetP(float y){Percent=y/100.0;};
	virtual ~ForceConst();
	friend std::ofstream & operator<<(std::ofstream &, ForceConst & );
	friend std::ifstream & operator>>(std::ifstream &, ForceConst & );
};

#endif /* FORCECONST_H_ */
