/*
 * Polarization.h
 *
 *  Created on: Dec 27, 2011
 *      Author: marchi
 */

#ifndef POLARIZATION_H_
#define POLARIZATION_H_

#include "Gridn.h"
#include "Rho.h"
#include "MyDipole.h"
#include "MyUtilClass.h"

using namespace DVECT;

namespace PolarizationNS {
inline float MyTriDiff(float g, float cubi){return((g-cubi>=0.0)? cubi-g+1 : g-cubi+1);};
class Polarization;
void MyDensity(Rho &,Polarization &,MyDipole &);

inline vector<float> MyTriLinear(Real gx, Real gy, Real gz){
	vector<float> dd(VERTEX);
	const rvec * cube=&Gridn<DIM>::Getcube()[0];
	for(int n=0;n<VERTEX;n++)
		dd[n]=MyTriDiff(gx,cube[n][XX])*MyTriDiff(gy,cube[n][YY])*MyTriDiff(gz,cube[n][ZZ]);
	return dd;
}


class Polarization: public Gridn<DIM> {
	int times;
	Dvect Mdip;
public:
	Polarization():Gridn<DIM>::Gridn(){
		std::string gh("Polar-output.bin");
		name=gh;
	};
	Dvect GetTotDip(){return Mdip;}
	Dvect & SetDip(){return Mdip;}
	Dvect Integrate() const{
		Dvect P_t=0.0;
		for(int i=0;i< static_cast<int>(nnx);i++)
			for(int j=0;j<static_cast<int>(nny);j++)
				for(int k=0;k<static_cast<int>(nnz);k++)
					for(int o=0;o<DIM;o++)
						P_t[o]+=(*this)[i][j][k][o];
		for(int o=0;o<DIM;o++) P_t[o]*=(this->getDV());
		return P_t;
	}
	void Density(MyDipole &);
	Polarization & operator=(const double);
	Polarization & operator=(const Polarization & y){
		Gridn<DIM>::operator=(y);
		return *this;
	}

	Polarization & operator+=(Polarization & y){
		Gridn<DIM>::operator+=(y);
		Mdip+=y.Mdip;
		times++;
		return *this;
	}
	Polarization & operator*=(const double y){
		Gridn<DIM>::operator*=(y);
		Mdip*=y;
		times++;
		return *this;
	}
	Polarization getavg(){
		try{
			if(!times) throw "Cannot take average if count is zero ";
		}
		catch (const char * s){
			std::cout << s << std::endl;
			exit(1);
		}

		Polarization temp(*this);
		temp/=static_cast<double>(times);
		return temp;
	}
	void getavgthis(){
		try{
			if(!times) throw "Cannot take average if count is zero ";
		}
		catch (const char * s){
			std::cout << s << std::endl;
			exit(1);
		}
		*this/=static_cast<double>(times);
	}

	virtual ~Polarization();
};

} /* namespace PolarizationNS */
#endif /* POLARIZATION_H_ */
