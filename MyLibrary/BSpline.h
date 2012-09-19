/*
 * BSpline.h
 *
 *  Created on: Jul 21, 2011
 *      Author: marchi
 */

#ifndef BSPLINE_H_
#define BSPLINE_H_
#include <vector>
using std::vector;

class BSpline {
protected:
	const int nders;
	static int order;
	vector<double> * spline;
	void Init(const double & w){
		spline[0][order-1]=0;
		spline[0][0]=1.0-w;
		spline[0][1]=w;
	};
	void OnePass(const double & w ,const int k){
		int o=k-1;
		double div;
		vector<double>::iterator c=spline[0].begin();
		div=1.0/ double (k-1);
		c[k-1]=c[k-2]*div*w;

	    for(int j=1;j<=k-2;j++){
	    	c[k-j-1]=div*((w+j)*c[k-j-2]+(k-j-w)*c[k-j-1]);
	    }
	    c[0]*=div*(1-w);
	}
	void Diff(){
		vector<double>::iterator c=spline[0].begin();
		vector<double>::iterator d=spline[1].begin();
		int j;
		d[0]=-c[0];
		for(int j=2;j<=order;j++){
			d[j-1]=c[j-2]-c[j-1];
		}
	}
public:
	BSpline();
	void Fill(const double & w);
	static void SetOrder(const int & o){order=o;};
	vector<double> & getSpline(){return spline[0];}
	vector<double> & getDSpline(){return spline[1];}
	virtual ~BSpline();
};

#endif /* BSPLINE_H_ */
