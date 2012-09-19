/*
 * EpsilonVor.h
 *
 *  Created on: Mar 29, 2012
 *      Author: marchi
 */

#ifndef EPSILONVOR_H_
#define EPSILONVOR_H_
#include "Epsilon.h"
#include <sstream>
#include <numeric>

namespace EpsilonNS {

class EpsilonVor: public Epsilon {
protected:
	static int nresid;
	vector<Matrix> D;
	vector<Matrix> G;
	vector<Matrix> D0;
	vector<Matrix> G0;
	vector<Matrix> eps;
	vector<Dvect> E;
	vector<Dvect> M;
	vector<Dvect> M0;
	vector<Dvect> E0;
	vector<string> ResLabel;
	vector<double> VoroVol;
	virtual void WriteIt(std::ofstream & ofout) const;
	virtual void ReadIt(std::ifstream & fin);
public:
	EpsilonVor(){
		try{
			nresid=Residue::GetNNR();
			if(!nresid) throw " Number of Residues not found by EpsilonVor. Allocation differred  ";
		}
		catch(const char * s){
			std::cout << s << std::endl;
			return;
		}
		Allocate();
		for(int o=0;o<nresid;o++) ResLabel[o]=Residue::GetLabel(o);
	};
	virtual EpsilonVor & operator-=(const Epsilon &);
	virtual EpsilonVor & operator=(const Epsilon &);


	virtual void TestConvergence(int);
	virtual Matrix TestD(int);
	virtual Matrix TestG(int);
	virtual Dvect TestM(int);
	virtual Dvect TestE(int);
	virtual ~EpsilonVor(){
		Deallocate();
	};
	virtual int GetNresid(){return nresid;};

	virtual bool operator()(Polarization &, Field &, Matrix &, Voronoi *);
	void testDipole(int,Dvect);

	virtual void Allocate(){
		Matrix zeroM=0.0;
		Dvect zeroV=0.0;
		D=vector<Matrix>(nresid,zeroM);
		G=vector<Matrix>(nresid,zeroM);
		eps=vector<Matrix>(nresid,zeroM);
		E=vector<Dvect>(nresid,zeroV);
		M=vector<Dvect>(nresid,zeroV);
		D0=vector<Matrix>(nresid,zeroM);
		G0=vector<Matrix>(nresid,zeroM);
		E0=vector<Dvect>(nresid,zeroV);
		M0=vector<Dvect>(nresid,zeroV);
		VoroVol=vector<double>(nresid,0.0);
		ResLabel=vector<string>(nresid);
	}
	virtual void Deallocate(){
		D.clear();
		G.clear();
		D0.clear();
		G0.clear();
		E.clear();
		M.clear();
		E0.clear();
		M0.clear();
		ResLabel.clear();
	}
};

} /* namespace EpsilonNS */
#endif /* EPSILONVOR_H_ */
