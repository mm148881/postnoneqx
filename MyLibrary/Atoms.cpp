/*
 * Atoms.cpp
 *
 *  Created on: May 20, 2011
 *      Author: marchi
 */
#include "Atoms.h"

Atoms::Atoms(const int nind) {
	nr=nind;
	if(nr) {
		x=new rvec [nr];
		xa=new rvec [nr];
	}
	else {x=NULL;xa=NULL;}
}
Atoms::Atoms(const AtomIndex & id) {
	nr=id.getN();
	if(nr) x=new rvec [nr];

	else x=NULL;
}
Atoms::Atoms(const Atoms & y) {
	nr=y.nr;
	x=new rvec [y.nr];
	for(int i=0;i<nr;i++) {
		x[i][0]=y.x[i][0];
		x[i][1]=y.x[i][1];
		x[i][2]=y.x[i][2];
	}
	Mt(y.Mt);
}

Atoms::~Atoms() {
	if(x) delete [] x;
	if(xa) delete [] xa;
	nr=0;
}
Atoms & Atoms::operator()(const int natoms){
	nr=natoms;
	if(!x) delete [] x;
	if(!xa) delete [] xa;
	x=new rvec [nr];
	xa=new rvec [nr];
	return *this;
};

Atoms & Atoms::operator=(const Atoms & y){
	try{
		if(nr != y.nr) throw "Cannot copy Atoms of different size ";
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
	return *this;
}
Atoms & Atoms::operator=(const real a){
	for(int i=0;i<nr;i++) {
		x[i][0]=a;
		x[i][1]=a;
		x[i][2]=a;
	}
	return *this;
}
void Atoms::Rot(const matrix coR){
	rvec t;
	for(int i=0;i<nr;i++){
		for(int o=0;o<DIM;o++) t[o]=coR[o][0]*x[i][0]+coR[o][1]*x[i][1]+coR[o][2]*x[i][2];
		for(int o=0;o<DIM;o++) x[i][o]=t[o];
	}
}
void Atoms::Tra(const rvec y){
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++) x[i][o]+=y[o];
}
void Atoms::setCoord(const Metric & Mt_in, const rvec * x0, const AtomIndex & id){
	try{
		if(nr != id.getN()) throw "wrong number of coordinates: Class Atoms";
	}
	catch(const char * s){
		std::cout<< s << std::endl;
		exit(1);
	}
	Mt(Mt_in);
	for(int i=0;i<id.getN();i++)
		for(int j=0;j<DIM;j++){
			x[i][j]=x0[id[i]][j];
		}
}
void Atoms::setCoord(const Metric & Mt_in, const rvec * x0){
	Mt(Mt_in);
	for(int i=0;i<nr;i++)
		for(int j=0;j<DIM;j++)
			x[i][j]=x0[i][j];
}
void Atoms::setMT(const Metric & Mt_in){
		Mt(Mt_in);
	}

void Atoms::doCOtoOC(){
	try{
	if(!x) throw "Atoms class atom coordinates undefined";
	}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	if(!xa) xa=new rvec [nr];
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++) xa[i][o]=Mt.getOC()[o][0]*x[i][0]+Mt.getOC()[o][1]*x[i][1]+
			Mt.getOC()[o][2]*x[i][2];
}
Atoms Atoms::COtoOC(){
	try{
	if(!x) throw "Atoms class atom coordinates undefined";
	}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	Atoms other(nr);
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++) other.x[i][o]=Mt.getOC()[o][0]*x[i][0]+Mt.getOC()[o][1]*x[i][1]+Mt.getOC()[o][2]*x[i][2];
	return other;
}
Atoms Atoms::OCtoCO(){
	try{
	if(!x) throw "Atoms class atom coordinates undefined";
	}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	Atoms other(nr);
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++) other.x[i][o]=Mt.getCO()[o][0]*x[i][0]+Mt.getCO()[o][1]*x[i][1]+Mt.getCO()[o][2]*x[i][2];
	return other;
}
Atoms Atoms::Shift(const rvec dx){
	Atoms tmp1(nr);
	tmp1=COtoOC();
	tmp1.Tra(dx);
	return tmp1;
}
real Atoms::Dist(const int i, const int j){
	try{
		if(!xa) throw "Atoms class atom reciprocal coordinates undefined";
	}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	rvec xc,xd;
	const matrix & co=Mt.getCO();
	for(int o=0;o<DIM;o++) {
		xd[o]=xa[i][o]-xa[j][o];
		xd[o]=xd[o]-rint(xd[o]);
	}
	for(int o=0;o<DIM;o++) xc[o]=co[o][0]*xd[0]+co[o][1]*xd[1]+co[o][2]*xd[2];

	return sqrt(xc[XX]*xc[XX]+xc[YY]*xc[YY]+xc[ZZ]*xc[ZZ]);
}
Atoms::plane Atoms::VectorDist(const int i, const int j){
	try{
		if(!xa) throw "Atoms class atom reciprocal coordinates undefined";
	}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	rvec xc,xd;
	const matrix & co=Mt.getCO();
	for(int o=0;o<DIM;o++) {
		xd[o]=xa[i][o]-xa[j][o];
		xd[o]=xd[o]-rint(xd[o]);
	}
	real frac=0.5;
	for(int o=0;o<DIM;o++) xc[o]=frac*(co[o][0]*xd[0]+co[o][1]*xd[1]+co[o][2]*xd[2]);
	ax=xc;
	ax=xc[XX]*xc[XX]+xc[YY]*xc[YY]+xc[ZZ]*xc[ZZ];
	ax=j;
	return ax;
}
void Atoms::pbc(){
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++)
			xa[i][o]=xa[i][o]-rint(xa[i][o]);
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++) x[i][o]=Mt.getCO()[o][0]*xa[i][0]+Mt.getCO()[o][1]*xa[i][1]+
			Mt.getCO()[o][2]*xa[i][2];

}
