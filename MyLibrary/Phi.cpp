/*
 * Phi.cpp
 *
 *  Created on: May 24, 2011
 *      Author: marchi
 */
#include "Phi.h"
void Phi::GetPhi(const Rho & myrho){

	unsigned int nzp=nnz/2+1;
	size_t align=sizeof(Complex);

	array3<double> fi(nnx,nny,nnz,align);
	array3<Complex> gi(nnx,nny,nzp,align);

	rcfft3d Forward(nnz,fi,gi);
	crfft3d Backward(nnz,gi,fi);

	int nfx,nfy,nfz;
	nfx=(nnx % 2 == 0)? nnx/2: nnx/2+1;
	nfy=(nny % 2 == 0)? nny/2: nny/2+1;
	nfz=(nnz % 2 == 0)? nnz/2: nnz/2+1;

	fi=myrho;

	Forward.fft(fi,gi);
	double mw1,mw2,mw3,mw,fact;
	int ia,ja,ka;
	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);

	for(int i=0;i<nx0;i++){
		ia=(i<nfx)?i : i-nx0;
		for(int j=0;j<ny0;j++){
			ja=(j<nfy)?j : j-ny0;
			for(int k=0;k<nz0/2+1;k++){
				ka=(k<nfz)?k : k-nz0;
				mw1=oc[XX][XX]*ia+oc[XX][YY]*ja+oc[XX][ZZ]*ka;
				mw2=oc[YY][XX]*ia+oc[YY][YY]*ja+oc[YY][ZZ]*ka;
				mw3=oc[ZZ][XX]*ia+oc[ZZ][YY]*ja+oc[ZZ][ZZ]*ka;
				mw=mw1*mw1+mw2*mw2+mw3*mw3;
				if(!j && !j && !k) fact=0.0;
				else fact=1.0/(mw*M_PI);
				gi[i][j][k]=fact*gi[i][j][k];
			}
		}
	}
	Backward.fftNormalized(gi,fi);
	*this=fi;
}
Phi::Phi(const Rho & myrho){
	GetPhi(myrho);
}
Phi & Phi::operator=(const Grid<1> & x){
	Grid<1> & My=*this;
	My=x;
	return *this;
}

Phi & Phi::operator=(array3<Complex> & x){
	Grid<1> & My=*this;
	My=x;
	return *this;
}
Phi & Phi::operator=(array3<double> & x){
	Grid<1> & My=*this;
	My[0]=x;
	return *this;
}

Phi & Phi::operator=(const double &  x){
	Grid<1> & My=*this;
	My=x;
	return *this;
}

Phi & Phi::operator()(const Rho & myrho){
	GetPhi(myrho);
	return *this;
}
