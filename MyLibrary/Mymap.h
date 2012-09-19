/*
 * Mymap.h
 *
 *  Created on: Oct 29, 2011
 *      Author: second
 */

#ifndef MYMAP_H_
#define MYMAP_H_

#include <map>

struct Tier{
private:
	static int mx,my,mz;
public:
	Tier(){};
	Tier(int cx,int cy,int cz){mx=cx;my=cy;mz=cz;};
	int nx;
	int ny;
	int nz;
	Tier & operator()(int cx,int cy,int cz){nx=cx;ny=cy;nz=cz;return *this;};
	bool operator<(const Tier & a) const{
		int my_int=mz*(my*nx+ny)+nz;
		int a_int=mz*(my*a.nx+a.ny)+a.nz;
		return my_int < a_int;
	}
	bool operator>(const Tier & a) const{
		return nx>a.nx && ny > a.ny && nz > a.nz;
	}
	bool operator==(const Tier & a) const{
		return nx == a.nx && ny == a.ny && nz == a.nz;
	}
};

template <class T1, class T2>
class Mymap: public std::map<T1,T2> {
public:
	Mymap();
	void insert(std::pair<T1,T2>);
	virtual ~Mymap();
};

#include "Mymap.cpp"
#endif /* MYMAP_H_ */
