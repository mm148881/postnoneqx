/*
 * Field.cpp
 *
 *  Created on: Dec 22, 2011
 *      Author: marchi
 */

#include "Field.h"

#define ORDER 4
array3<double> * Field::ro0=NULL;
array3<Complex> * Field::rok0=NULL;
array3<Complex> * Field::dgi=NULL;
Field & Field::operator+=(Field & y){
	Gridn<DIM>::operator+=(y);
	times++;
	return *this;
}

Field & Field::operator()(Grid<1> & y){
	array4<double> dfi0(DIM,nnx,nny,nnz);
	times=1;
	int nfx,nfy,nfz;
	nfx=(nnx % 2 == 0)? nnx/2: nnx/2+1;
	nfy=(nny % 2 == 0)? nny/2: nny/2+1;
	nfz=(nnz % 2 == 0)? nnz/2: nnz/2+1;
	unsigned int nzp=nnz/2+1;
	size_t align=sizeof(Complex);
	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);

	array3<double> & ro=*ro0;
	array3<Complex> & rok=*rok0;

	rcfft3d Forward(nnz,ro,rok);
	crfft3d Backward(nnz,rok,ro);

	for(int i=0;i<nx0;i++)
		for(int j=0;j<ny0;j++)
			for(int k=0;k<nz0;k++){
				ro[i][j][k]=y[0][i][j][k];
			}
	Forward.fft(ro,rok);

	if(filter) {
		array3<Complex> & Mfilter=*filterrc;

		for(int i=0;i<nx0;i++)
			for(int j=0;j<ny0;j++)
				for(int k=0;k<nz0/2+1;k++){
					rok[i][j][k]=rok[i][j][k]*Mfilter[i][j][k];
				}
		}
	double mw1,mw2,mw3,mw,fact,fact1;
	int ia,ja,ka;

	Complex imag(0.0,1.0), zero(0,0);
	double mysinx, mysiny,mysinz;
	DiffCoeffs<ORDER> My(nx0,ny0,nz0,co[0][0],co[1][1],co[2][2]);

	for(int i=0;i<nx0;i++){
		ia=(i<nfx)?i : i-nx0;
		mysinx=My.coeffs(XX,ia);
		for(int j=0;j<ny0;j++){
			ja=(j<nfy)?j : j-ny0;
			mysiny=My.coeffs(YY,ja);
			for(int k=0;k<nz0/2+1;k++){
				ka=(k<nfz)?k : k-nz0;
				mysinz=My.coeffs(ZZ,ka);

				mw1=oc[XX][XX]*ia+oc[XX][YY]*ja+oc[XX][ZZ]*ka;
				mw2=oc[YY][XX]*ia+oc[YY][YY]*ja+oc[YY][ZZ]*ka;
				mw3=oc[ZZ][XX]*ia+oc[ZZ][YY]*ja+oc[ZZ][ZZ]*ka;
				mw=mw1*mw1+mw2*mw2+mw3*mw3;
				if(!j && !j && !k) fact=0.0;
				else fact=1.0/(mw*M_PI);
				dgi[XX][i][j][k]=mysinx*fact*imag*rok[i][j][k];
				dgi[YY][i][j][k]=mysiny*fact*imag*rok[i][j][k];
				dgi[ZZ][i][j][k]=mysinz*fact*imag*rok[i][j][k];
			}
		}
	}

	Backward.fftNormalized(dgi[XX],dfi0[XX]);
	Backward.fftNormalized(dgi[YY],dfi0[YY]);
	Backward.fftNormalized(dgi[ZZ],dfi0[ZZ]);
	for(unsigned int i=0;i<nnx;i++)
		for(unsigned int j=0;j<nny;j++)
			for(unsigned int k=0;k<nnz;k++){
				(*this)[i][j][k][XX] = -dfi0[XX][i][j][k];
				(*this)[i][j][k][YY] = -dfi0[YY][i][j][k];
				(*this)[i][j][k][ZZ] = -dfi0[ZZ][i][j][k];
			}

	return *this;
}

Field::~Field() {
	// TODO Auto-generated destructor stub
}

