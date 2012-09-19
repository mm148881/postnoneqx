/*
 * Metric.cpp
 *
 *  Created on: May 25, 2011
 *      Author: marchi
 */

#include "Metric.h"
const matrix Metric::Idtity={{1,0,0},{0,1,0},{0,0,1}};
Metric::Metric(){
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) {
			co[i][j]=Idtity[i][j];
			oc[i][j]=Idtity[i][j];
		}

}
Metric::Metric(const matrix & co_in) {
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) co[i][j]=co_in[i][j];
	invertCO();
}
Metric::Metric(const Metric & Mt_in) {
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) {
			co[i][j]=Mt_in.co[i][j];
			oc[i][j]=Mt_in.oc[i][j];
		}

}
Metric Metric::operator/(const double a){
	Metric tmp=*this;
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) {
			tmp.co[i][j]=tmp.co[i][j]/(a);
			tmp.oc[i][j]=tmp.oc[i][j]/(a);
		}
	return tmp;
}
Metric & Metric::operator()(const Metric & Mt_in) {
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) {
			co[i][j]=Mt_in.co[i][j];
			oc[i][j]=Mt_in.oc[i][j];
		}
	return *this;
}
Metric & Metric::operator=(const Metric & Mt_in) {
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) {
			co[i][j]=Mt_in.co[i][j];
			oc[i][j]=Mt_in.oc[i][j];
		}
	return *this;
}
Metric & Metric::operator+=(const Metric & Mt_in) {
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) {
			co[i][j]+=Mt_in.co[i][j];
			oc[i][j]+=Mt_in.oc[i][j];
		}
	return *this;
}
Metric & Metric::operator()(const matrix & co_in){
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) co[i][j]=co_in[i][j];
	invertCO();
	return *this;
}

Metric::~Metric() {

}
void Metric::invertCO()
{
  /* Save some time by assuming lower right part is zero */
  real tmp=1.0/(co[XX][XX]*co[YY][YY]*co[ZZ][ZZ]);
  oc[XX][XX]=co[YY][YY]*co[ZZ][ZZ]*tmp;
  oc[YY][XX]=0;
  oc[ZZ][XX]=0;
  oc[XX][YY]=-co[XX][YY]*co[ZZ][ZZ]*tmp;
  oc[YY][YY]=co[XX][XX]*co[ZZ][ZZ]*tmp;
  oc[ZZ][YY]=0;
  oc[XX][ZZ]=(co[XX][YY]*co[YY][ZZ]-co[YY][YY]*co[XX][ZZ])*tmp;
  oc[YY][ZZ]=-co[YY][ZZ]*co[XX][XX]*tmp;
  oc[ZZ][ZZ]=co[XX][XX]*co[YY][YY]*tmp;
}
