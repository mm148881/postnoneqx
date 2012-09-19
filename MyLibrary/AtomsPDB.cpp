/*
 * AtomsPDB.cpp
 *
 *  Created on: Oct 26, 2011
 *      Author: marchi
 */

#include "AtomsPDB.h"
int AtomsPDB::count=0;

AtomsPDB::AtomsPDB(t_atoms * atoms,AtomIndex & id):Atoms::Atoms(id){
	fout=NULL;
	MyId=new AtomIndex;
	*MyId=id;
    init_t_atoms(&useatoms,atoms->nr,FALSE);
    delete [] useatoms.resinfo;
    useatoms.resinfo = atoms->resinfo;
    for(int i=0;(i<id.getN());i++) {
        useatoms.atomname[i]=atoms->atomname[id[i]];
        useatoms.atom[i]=atoms->atom[id[i]];
        useatoms.nres=max(useatoms.nres,useatoms.atom[i].resind+1);
    }
    useatoms.nr=id.getN();
}

void AtomsPDB::setCoord(const Metric & Mt_in, const rvec * x0){
	Mt(Mt_in);
	for(int i=0;i<MyId->getN();i++)
		for(int j=0;j<DIM;j++){
			x[i][j]=x0[(*MyId)[i]][j];
		}
}
AtomsPDB & AtomsPDB::operator=(const real y){
	for(int i=0;i<nr;i++) {
		x[i][0]=y;
		x[i][1]=y;
		x[i][2]=y;
	}
	return *this;
}
AtomsPDB & AtomsPDB::operator=(const AtomsPDB & y){
	try{
		if(nr != y.nr) throw "Cannot copy AtomsPDB of different size ";
	}
	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	for(int i=0;i<nr;i++) {
		x[i][0]=y.x[i][0];
		x[i][1]=y.x[i][1];
		x[i][2]=y.x[i][2];
	}
	Mt(y.Mt);

	if(!MyId) MyId=new AtomIndex;
	else {
		delete MyId;
		MyId=new AtomIndex;
	}
	*MyId=*(y.MyId);
    init_t_atoms(&useatoms,y.useatoms.nr,FALSE);
    delete [] useatoms.resinfo;
    useatoms.resinfo = y.useatoms.resinfo;
    for(int i=0;i<y.useatoms.nr;i++) {
        useatoms.atomname[i]=y.useatoms.atomname[i];
        useatoms.atom[i]=y.useatoms.atom[i];
        useatoms.nres=max(useatoms.nres,y.useatoms.atom[i].resind+1);
    }
    useatoms.nr=y.useatoms.nr;
	return *this;
}
AtomsPDB &  AtomsPDB::operator+=(const Atoms & y){
	try{
		if(nr != y.getNR()) throw "Cannot accumulate AtomsPDB of different size ";
	}
	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	Mt+=y.getMt();
	for(int i=0;i<nr;i++) {
		x[i][0]+=y[i][0];
		x[i][1]+=y[i][1];
		x[i][2]+=y[i][2];
	}
	count++;
	Metric ui=y.getMt();
	return *this;
}
void AtomsPDB::Print(const char * title,int grps){
	rvec * x0=new rvec [nr];
	real cnt=1.0;
	if(count) cnt=static_cast<real> (count);
	Metric tmp=Mt/static_cast<double> (cnt);
	for(int i=0;i<nr;i++) for(int o=0;o<DIM;o++) x0[i][o]=x[i][o]/cnt;
	matrix & co=const_cast<matrix &> (tmp.getCO());
	write_pdbfile(fout,title,&useatoms,x0,0,co,' ',grps,NULL,TRUE);
}

AtomsPDB::~AtomsPDB() {
	// TODO Auto-generated destructor stub
}

