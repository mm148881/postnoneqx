/*
 * Rho.h
 *
 *  Created on: May 24, 2011
 *      Author: marchi
 */

#ifndef RHO_H_
#define RHO_H_


#include <string>
#include "typedefs.h"
#include "physics.h"
#include "Charges.h"
#include "Grid.h"
#include <typeinfo>
#include <cstring>
#include <set>
using std::multiset;
using std::set;
using namespace Array;


struct DataDensity{
	int n;
	double value;
};
struct DensityComp {
  bool operator() (const DataDensity & lhs, const DataDensity & rhs) const
  {return lhs.n<rhs.n;}
};

class Rho : public Grid<1> {
protected:
public:
	Rho():Grid<1>::Grid(){std::string gh("RhoOutput.bin");name=gh;};
	Rho & operator=(const Grid<1> &);
	Rho & operator=(const double &);
	void Density(Charges & );
	void Density1(Charges & );
	void Symmetrify(rvec & x);
	double Sum(){
		double tot=0.0;
		for(unsigned int i=0;i<this->Size();i++) tot+=this->v[i];
		return tot;
	};
	void Linear(Rho,Rho);
	virtual void Accumulate(){};
	virtual ~Rho(){};
	friend std::ofstream & operator<<(std::ofstream &, Rho & );
	friend std::ifstream & operator>>(std::ifstream &, Rho & );
};

#endif /* RHO_H_ */
