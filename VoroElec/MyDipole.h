/*
 * MyDipole.h
 *
 *  Created on: Dec 27, 2011
 *      Author: marchi
 */

#ifndef MYDIPOLE_H_
#define MYDIPOLE_H_

#include <iostream>
#include <Charges.h>
#include <vector>
#include "MyUtilClass.h"
using std::vector;
using DVECT::Dvect;

typedef double dvect[DIM];
typedef::real Myreal;

struct Dip{
	dvect d; // The dipole moment
	dvect x; // the center of mass in Cartesian coordinates
	dvect xa; /// the center of mass in reduced coordinates
	Dip & operator=(const Dip & y){
		for(int i=0;i<DIM;i++) d[i]=y.d[i];
		for(int i=0;i<DIM;i++) x[i]=y.x[i];
		for(int i=0;i<DIM;i++) xa[i]=y.xa[i];
		return *this;
	}
		Dip & operator=(const double y){
		for(int i=0;i<DIM;i++) d[i]=y;
		for(int i=0;i<DIM;i++) x[i]=y;
		for(int i=0;i<DIM;i++) xa[i]=y;
		return *this;
	}
};
struct DipG{
	double q, tmass;
	vector<double> mass;
	vector<int> cidx;
	DipG & operator=(const DipG & y){
		mass=y.mass;
		cidx=y.cidx;
		q=y.q;
		tmass=y.tmass;
		return *this;
	}
};
class MyDipole: public Charges {
	int ndipoles;
	int natoms;
	Dip * dip;
	DipG * Conf;
	static vector<int> resid;
	static vector<int> residDip;
	static int nresid;
public:
	MyDipole():Charges(), ndipoles(0),natoms(0),dip(NULL),Conf(NULL){};
	MyDipole(const int,const int, const int *, const Myreal *, const Myreal *, const vector<int> *,int );
	MyDipole(const MyDipole &);
	void setCoord(const Metric &, const rvec *);
	virtual ~MyDipole();
	dvect & getDip(int i){return dip[i].d;};
	dvect & getX(int i){return dip[i].x;};
	dvect & getXa(int i){return dip[i].xa;};
	int getNdip(){return ndipoles;};
	int getNatom(){return natoms;};
	bool test(){return (dip)?TRUE:FALSE;}
};

#endif /* MYDIPOLE_H_ */
