/*
 * Residue.h
 *
 *  Created on: Jan 18, 2012
 *      Author: marchi
 */

#ifndef RESIDUE_H_
#define RESIDUE_H_
#include <vector>
#include <string>
#include <iostream>
#include "typedefs.h"
#include "AtomIndex.h"


using std::cout;
using std::endl;

using std::vector;
using std::string;
class Residue {
protected:
	static int nnr;
	static int nr;
	static vector<int> * ind;
	static vector<int> * ind0;
	static int * resid;
	static int * resid0;
	static string * label;
public:
	Residue();
	static void SetUp(const int, const int, const vector<int> *,const vector<int> *, const string *);
	static void Select(const AtomIndex & cidx);
	static int GetNR(){return nr;};
	static int GetNNR(){return nnr;};
	static int ResiduePid(const int n){
		try{
			if(n > nr || n <0) throw " ResiduePid Index inconsistent with number of atoms ";
		}
		catch (const char * s) {
			cout << s << endl;
			exit(1);
		}
		return resid[n];
	};
	static int ResiduePid0(const int n){
		try{
			if(n > nr || n <0) throw " ResiduePid Index inconsistent with number of atoms ";
		}
		catch (const char * s) {
			cout << s << endl;
			exit(1);
		}
		return resid0[n];
	};
	static int Pid(const int n){
		try{if(n > nnr || n <0) throw "\n Pid Index inconsistent with number of atoms ";}
		catch (const char * s) {cout << s << endl;exit(1);}
		int out=ind[n].size()?resid[ind[n][0]]:-1;
		return out;
	}
	static vector<int> & AtomsPids(const int n){
		try{if(n > nnr || n <0) throw "\nAtomPid Index inconsistent with number of atoms ";}
		catch (const char * s) {cout << s << endl;exit(1);}
		return ind[n];
	};
	static string GetLabel(int n){
		try{
			if(n >= nnr) throw "\n Asking for label of a residue number larger than the total ";
		} catch(const char * s){
			std::cout << s << std::endl;
			exit(-1);
		}
		return label[n];
	}
	virtual ~Residue();
};

#endif /* RESIDUE_H_ */
