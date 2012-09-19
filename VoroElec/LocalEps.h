/*
 * LocalEps.h
 *
 *  Created on: Jan 18, 2012
 *      Author: marchi
 */

#ifndef LOCALEPS_H_
#define LOCALEPS_H_
#include "Residue.h"
#include <string>
#include <vector>
//#include "maths.h"
#include "MyUtilClass.h"
#include <iomanip>
#include "stdafx.h"
#include <math.h>
#include "linalg.h"
#include <fstream>

using namespace std;
using namespace alglib;
using std::string;
using std::vector;
using namespace DVECT;
using namespace MATRIX;

class LocalEps: public Residue {
	Dvect * P;
	Dvect * E;

	Dvect M_avg;
	Dvect * P_avg;
	Dvect * E_avg;
	Matrix * PM_avg;
	Matrix * EM_avg;
	Matrix * eps;
	Matrix * D_fl;
	Matrix * G_fl;
	Matrix * Gm1_fl;
	Dvect P0_avg,E0_avg;
	Matrix PM0_avg,EM0_avg,eps0,D0_fl,G0_fl,Gm10_fl;

	int * count;
	int TotalCount;
public:
	LocalEps(): Residue(), P(NULL), E(NULL), eps(NULL), P_avg(NULL), E_avg(NULL),
	PM_avg(NULL), EM_avg(NULL), D_fl(NULL), G_fl(NULL), Gm1_fl(NULL), count(NULL), TotalCount(0) {};
	LocalEps(LocalEps & y): Residue(), P(NULL), E(NULL), eps(NULL), P_avg(NULL), E_avg(NULL),
	PM_avg(NULL), EM_avg(NULL), count(NULL), D_fl(NULL), G_fl(NULL), Gm1_fl(NULL), TotalCount(0) {copy(y);};
	LocalEps & operator=(LocalEps & y){copy(y);return *this;};
	LocalEps & operator()(LocalEps & y){copy(y);return *this;};

	void copy(LocalEps & y ){
		SetUp();
		P=y.P;
		E=y.E;
		P_avg=y.P_avg;
		E_avg=y.E_avg;
		PM_avg=y.PM_avg;
		EM_avg=y.EM_avg;
		eps=y.eps;
		count=y.count;
		M_avg=y.M_avg;
		TotalCount=y.TotalCount;
		G_fl=y.G_fl;
		D_fl=y.D_fl;
		Gm1_fl=y.Gm1_fl;
	}
	double getP(int i,int j){return P[i][j];};
	void SetUp();
	void Set(const int i,const Dvect & x,const Dvect & y){
		try {if(i>=nr) throw "Goes above dimensions set";}
		catch(const char * s){std::cout << s << std::endl;exit(1);}
		P[i]=x;
		E[i]=y;
	};
	void Accumulate(Dvect & M){TotalCount++; M_avg+=M;};
	void Accumulate(const int i,const Dvect & p,const Dvect & e, Dvect & M){
		try {if(i>=nr) throw "Goes above dimensions set";}
		catch(const char * s){std::cout << s << std::endl;exit(1);}
		P_avg[i]+=p;
		E_avg[i]+=e;

		PM_avg[i]+=p%M;
		EM_avg[i]+=e%M;

/*
		PM_avg[i]+=p%p;
		EM_avg[i]+=e%p;
*/
		count[i]++;
	};
	void Avg();
	void Zero();
	void Print();
	Matrix & GetG_fl(int i) const {return G_fl[i];}

	virtual ~LocalEps(){
		if(P) delete [] P;
		if(E) delete [] E;
		if(P_avg) delete [] P_avg;
		if(E_avg) delete [] E_avg;
		if(PM_avg) delete [] PM_avg;
		if(EM_avg) delete [] EM_avg;
		if(G_fl) delete [] G_fl;
		if(D_fl) delete [] D_fl;
		if(eps) delete [] eps;
		if(count) delete [] count;
	}
	friend std::ostream & operator<<(std::ostream &, LocalEps & );
	friend std::ofstream & operator<<(std::ofstream &, LocalEps & );
};

#endif /* LOCALEPS_H_ */
