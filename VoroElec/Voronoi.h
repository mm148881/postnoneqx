/*
 * Voronoi.h
 *
 *  Created on: Jun 20, 2011
 *      Author: marchi
 */

#ifndef VORONOI_H_
#define VORONOI_H_
#include <string>
#include <vector>
#include <list>
#include "Atoms.h"
#include "ResidueCM.h"
#include <cstring>
#include "Array.h"
#include <fstream>
#include <iomanip>
#include "VoronoiPrint.h"

using namespace Array;

using namespace std;
#include "voro++.hh"
using namespace voro;

using namespace std;
const int NNN=26;
const char h[]="Hxxxxx";


class Voronoi {
	static int nresid,nr,nc;
	static float time;
	static int * types;
	string * TypesName;
	double * Vol;
	array1<double> Vols;
	array2<double> area;
	container_periodic * Mycon;
	particle_order * porder;
	static string * label;
	vector<int> cindex;
	vector<int> * Neighs;
	vector<double> * Surface;
	vector<int> RealResidue;
	void gather(vector<int> & y);
	Voronoi();
public:
	Voronoi(int, char ***,bool);
	void Start(float, Atoms &);
	void getData();
	virtual ~Voronoi(){
		delete [] TypesName;
		delete [] Vol;
		delete [] Neighs;
		delete [] Surface;
		delete Mycon;
		delete porder;
	};
	string & name(int i) const {return label[i];}
	void setTypes(int,const AtomIndex [],const string *);
	double & getVolR(int n){return Vols[n];};
	int VoronoiCellID(const double & xc,const double & yc,const double & zc){
		double rx,ry,rz;
		int iD;
		try{if(!Mycon->find_voronoi_cell(xc,yc,zc,rx,ry,rz,iD))
			throw "Error while finding Voronoi cell ";}
		catch(const char * s) {std::cout << s << std::endl;exit(1);}
		return cindex[iD];
	}
	int VoronoiCellID(const double & xc,const double & yc,const double & zc,
			double & rx, double & ry, double & rz){
		int iD;
		try{if(!Mycon->find_voronoi_cell(xc,yc,zc,rx,ry,rz,iD))
			throw "Error while finding Voronoi cell ";}
		catch(const char * s) {std::cout << s << std::endl;exit(1);}
		return cindex[iD];
	}
	int getTypes(int i) {return types[i];};
	int getTypesRes(int nc) {
		int First=ResidueCM::getind(nc)[0];
		return getTypes(First);
	}
	string & getTypesName(int i) const {return TypesName[types[i]];};
	friend std::ofstream & operator<<(std::ofstream &, Voronoi  & );
};
#endif /* VORONOI_H_ */
