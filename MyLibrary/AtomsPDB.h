/*
 * AtomsPDB.h
 *
 *  Created on: Oct 26, 2011
 *      Author: marchi
 */

#ifndef ATOMSPDB_H_
#define ATOMSPDB_H_
#include "Atoms.h"
#include "macros.h"
#include "pdbio.h"

class AtomsPDB: public Atoms {
protected:
	t_atoms useatoms;
	AtomIndex * MyId;
	static int count;
	int ePDB;
	FILE * fout;
public:
	AtomsPDB() {ePDB=0;fout=NULL;MyId=NULL;nr=0; x=NULL; xa=NULL;};
	AtomsPDB(t_atoms *,AtomIndex &);
	void setCoord(const Metric &, const rvec *);

	void Print(const char *,int);
	void openstream(FILE * y){
		fout=y;
	}
	AtomsPDB & operator+=(const Atoms &);
	AtomsPDB & operator=(const AtomsPDB & );
	AtomsPDB & operator=(const real );
	virtual ~AtomsPDB();
};

#endif /* ATOMSPDB_H_ */
