/*
 * ResidueCM.h
 *
 *  Created on: Sep 21, 2011
 *      Author: marchi
 */

#ifndef RESIDUECM_H_
#define RESIDUECM_H_
#include <vector>
#include "Atoms.h"
#include "typedefs.h"
#include <string>
#include "Array.h"
#include <sstream>
#include <fstream>
#include <iomanip>

using namespace Array;

using std::string;

using std::vector;
class ResidueCM:public Atoms {
private:
	static int nnr;
	static vector<int> * ind;
	static string * label;
	static vector<int> * iD;
	rvec * xcm;
	rvec * xcma;
	void Allocate(int);
public:
	ResidueCM(): Atoms::Atoms(), xcm(NULL), xcma(NULL){};
	ResidueCM(const int, const int, const vector<int> *, const string *);
	ResidueCM(const int natoms, const int n): Atoms::Atoms(natoms) {nnr=n; ind=new vector<int> [nnr];xcm=new rvec [nnr];}
	ResidueCM(const rvec *);
	virtual ~ResidueCM();
	virtual void Allocate();
	ResidueCM & operator()(const rvec *);
	ResidueCM & operator()(const int );
	ResidueCM & operator()(const ResidueCM &);
	ResidueCM & operator()(const int, const int, const vector<int> *, const string *);
	ResidueCM & operator=(const float);
	ResidueCM & operator=(const ResidueCM &);
	ResidueCM & operator+=(const ResidueCM &);
	ResidueCM operator+(const ResidueCM &);
	rvec & operator[](const int i)const {return xcm[i];};
	void doCOtoOC();
	void ResidueSelect(const int,const AtomIndex *);
	static vector<int> & getind(int i) {return ind[i];}
	static int   getiD(int i) {return (*iD)[i];};
	static int   getiDSize(){return iD->size();};
	void setlabels();
	static string & getlabels(int i) {return label[i];};
	static int Size(){return nnr;};
	void SetCM(const rvec *);
	void SetCM();
	friend std::ostream & operator<<(std::ostream &, const ResidueCM  & );
	virtual void pbc();
};

#endif /* RESIDUECM_H_ */
