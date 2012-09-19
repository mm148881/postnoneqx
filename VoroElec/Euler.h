/*
 * Euler.h
 *
 *  Created on: Feb 28, 2012
 *      Author: marchi
 */

#ifndef EULER_H_
#define EULER_H_
#include "MyUtilClass.h"
#include <cmath>
using namespace DVECT;
using namespace MATRIX;


class Euler {
	double Alpha,Beta,Gamma;
	inline void set(Dvect y0){
		Dvect z1;
		if(y0 != 0.0) z1=y0.normal();
		else z1=y0;

		Dvect x(1.0,0.0,0.0);
		Dvect y(0.0,1.0,0.0);
		Dvect z(0.0,0.0,1.0);
		Dvect cross=z^z1;
		if(cross == 0.0) {
			Alpha=0.0;
			Beta=0.0;
			Gamma=0.0;
			return;
		} else cross.normalize();
		double arg1=cross*x;
		double arg0=cross*y;

		double arg2=z*z1;
		Alpha=atan2(arg0,arg1);
		Beta=acos(arg2);
		Gamma=0.0;
	}
public:
	Euler(): Alpha(0.0),Beta(0.0), Gamma(0.0){};
	Euler(const Euler & y){Alpha=y.Alpha;Beta=y.Beta;Gamma=y.Gamma;};
	Euler(const Dvect & y){set(y);};
	void operator()(const Dvect & y){Euler tmp(y);*this=tmp;};
	Euler operator=(Euler & y){Alpha=y.Alpha;Beta=y.Beta;Gamma=y.Gamma;return *this;};
	double getAlpha(){return Alpha;};
	double getBeta(){return Beta;};
	double getGamma(){return Gamma;};

	Matrix getMatrix(){
		Matrix m;
		double & f=Alpha;
		double & t=Beta;
		m[XX][XX]=cos(f);
		m[XX][YY]=sin(f);
		m[XX][ZZ]=0.0;
		m[YY][XX]=-cos(t)*sin(f);
		m[YY][YY]=cos(t)*cos(f);
		m[YY][ZZ]=sin(t);
		m[ZZ][XX]=sin(t)*sin(f);
		m[ZZ][YY]=-sin(t)*cos(f);
		m[ZZ][ZZ]=cos(t);
		return m;
	}

	inline const Dvect Rotate(const Dvect & My){
		Matrix Rot=this->getMatrix();
		Dvect tmp=0.0;
		for(int o=0;o<DIM;o++){
			tmp[XX]+=Rot[XX][o]*My[o];
			tmp[YY]+=Rot[YY][o]*My[o];
			tmp[ZZ]+=Rot[ZZ][o]*My[o];
		}
		return tmp;
	}
	inline Matrix Rotate(const Matrix & My){
		Matrix My1=My;
		Matrix Rot=this->getMatrix();
		Matrix Rot_t=Rot.Transpose();
		Matrix tmp=My1*Rot_t;
		Matrix tmp1=Rot*tmp;
		return tmp1;
	}
	virtual ~Euler();
};

#endif /* EULER_H_ */
