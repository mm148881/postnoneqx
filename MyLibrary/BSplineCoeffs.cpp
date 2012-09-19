/*
 * BSplineCoeffs.cpp
 *
 *  Created on: Jul 21, 2011
 *      Author: marchi
 */

#include "BSplineCoeffs.h"

void BSplineCoeffs::MakeCoeffs(Charges & q){
	try {
		if(!(nnx || nny || nnz)) throw "Cannot make cspline coefficients if grid dimension are not set ";
	}
	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	if(nr != q.getNR()){
		if(coeffs) delete [] coeffs;
		coeffs=new BSpline3D[nr];
		Atoms x=q.COtoOC();
		for(int i=0;i<nr;i++){
			double wX=double (nnx)*(x[i][XX]-rint(x[i][XX]));
			double wY=double (nny)*(x[i][YY]-rint(x[i][YY]));
			double wZ=double (nnz)*(x[i][ZZ]-rint(x[i][ZZ]));
			coeffs[i][XX].Fill(wX);
			coeffs[i][YY].Fill(wY);
			coeffs[i][ZZ].Fill(wZ);
		}
	}
}
BSplineCoeffs::BSplineCoeffs() {
	// TODO Auto-generated constructor stub

}

BSplineCoeffs::~BSplineCoeffs() {
	// TODO Auto-generated destructor stub
}

