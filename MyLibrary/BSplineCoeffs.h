/*
 * BSplineCoeffs.h
 *
 *  Created on: Jul 21, 2011
 *      Author: marchi
 */

#ifndef BSPLINECOEFFS_H_
#define BSPLINECOEFFS_H_

#include "Charges.h"
#include "BSpline.h"
#include "typedefs.h"
#include <cmath>

class BSplineCoeffs{
	typedef BSpline BSpline3D[DIM];
	BSpline3D * coeffs;
	int nnx,nny,nnz;
	int nr;
	BSplineCoeffs();
	void MakeCoeffs(Charges & x);
public:
	BSplineCoeffs(int nx0, int ny0, int nz0):
		nnx(nx0), nny(ny0), nnz(nz0), nr(0), coeffs(NULL){};
	void operator()(Charges & x,const int order){	BSpline::SetOrder(order);MakeCoeffs(x);};
	void operator()(Charges & x){MakeCoeffs(x);};
	virtual ~BSplineCoeffs();
};

#endif /* BSPLINECOEFFS_H_ */
