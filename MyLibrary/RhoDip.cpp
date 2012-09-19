/*
 * RhoDip.cpp
 *
 *  Created on: Jun 21, 2011
 *      Author: marchi
 */

#include "RhoDip.h"
#include <iostream>
#include <cmath>


RhoDip & RhoDip::operator=(const Grid<DIM> & x){
	Grid<DIM> & My=*this;
	My=x;
	return *this;
}
RhoDip & RhoDip::operator=(const double & x){
	Grid<DIM> & My=*this;
	My=x;
	return *this;
}
void RhoDip::Density(Dipoles & x){
	try{
		if(!x.test()) throw "Cannot compute density if Dipoles not initialized.";
	}
	catch(const char *s){
		std::cout << s << std::endl;
	}
	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);


	int natoms=x.getNR();
	array2<real> xa(natoms,DIM);

	real x1,y1,z1,r1,s1,t1,gx,gy,gz;
	int mx,my,mz;
	int ox,oy,oz,tx,ty,tz;
	double DV=RhoDip::getDV();
	if(strstr(typeid(*this).name(),"RhoDipAverage") == NULL) *this=0.0;
	else this->Accumulate();

	for (int i = 0; i < natoms; i++) {
		xa[i][XX] = oc[XX][XX] * x[i][XX] + oc[XX][YY] * x[i][YY] + oc[XX][ZZ] * x[i][ZZ];
		xa[i][YY] = oc[YY][XX] * x[i][XX] + oc[YY][YY] * x[i][YY] + oc[YY][ZZ] * x[i][ZZ];
		xa[i][ZZ] = oc[ZZ][XX] * x[i][XX] + oc[ZZ][YY] * x[i][YY] + oc[ZZ][ZZ] * x[i][ZZ];
	}
	for (int i = 0; i < natoms; i++) {
		x1 = xa[i][XX];
		y1 = xa[i][YY];
		z1 = xa[i][ZZ];
		r1 = static_cast<real> (nx0 * (x1 - rint(x1 - 0.5)));
		s1 = static_cast<real> (ny0 * (y1 - rint(y1 - 0.5)));
		t1 = static_cast<real> (nz0 * (z1 - rint(z1 - 0.5)));
		mx = static_cast<int>(r1);
		my = static_cast<int>(s1);
		mz = static_cast<int>(t1);
		gx = r1 - (real) mx;
		gy = s1 - (real) my;
		gz = t1 - (real) mz;
		rvec & chg=x.gdip(i);
		TriLinear(gx, gy, gz, chg);
		for (int n = 0; n < VERTEX; n++) {
			tx = static_cast<int>(cube[n][XX]);
			ty = static_cast<int>(cube[n][YY]);
			tz = static_cast<int>(cube[n][ZZ]);
			ox = (nx0 > mx + tx) ? mx + tx : mx + tx - nx0;
			oy = (ny0 > my + ty) ? my + ty : my + ty - ny0;
			oz = (nz0 > mz + tz) ? mz + tz : mz + tz - nz0;
			for(int o=0;o<DIM; o++)
				(*this)[o][ox][oy][oz] += static_cast<double>(qq[n][o])/DV;
		}
	}
	xa.Deallocate();
}
std::ofstream & operator<<(std::ofstream & fout, RhoDip & y){
	fout.write((char *) &y.nnx, sizeof y.nnx);
	fout.write((char *) &y.nny, sizeof y.nny);
	fout.write((char *) &y.nnz, sizeof y.nnz);
	fout.write((char *) &y.co, sizeof y.co);
	fout.write((char *) &y.oc, sizeof y.oc);
	fout.write((char *) &y[0][0][0][0], (sizeof y[0][0][0][0])*DIM*y.nnx*y.nny*y.nnz);
	return fout;
}
std::ifstream & operator>>(std::ifstream & fin, RhoDip & y){
	unsigned int nnx,nny,nnz;
	matrix co,oc;
	fin.read((char *) &nnx, sizeof nnx);
	fin.read((char *) &nny, sizeof nny);
	fin.read((char *) &nnz, sizeof nnz);
	fin.read((char *) &co, sizeof co);
	fin.read((char *) &oc, sizeof oc);
	y.set(co,nnx,nny,nnz);
	y.Allocate();
	fin.read((char *) &y[0][0][0][0], (sizeof y[0][0][0][0])*DIM*nnx*nny*nnz);
	return fin;
}
