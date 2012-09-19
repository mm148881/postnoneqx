/*
 * NeighListB.h
 *
 *  Created on: 9 Aug 2011
 *      Author: marchi
 */

#ifndef NEIGHLIST_H_
#define NEIGHLIST_H_
#include <list>

struct NonSortable{
	int n;
	NonSortable(const int & a){n=a;}
	NonSortable(){};
};
struct Sortable{
	int n;
	float dist;
	Sortable(){};
	Sortable(const float & a, const int & b){n=b;dist=a;};
	Sortable operator()(const float & a, const int & b){n=b;dist=a;return *this;};
	Sortable operator=(Sortable & x){n=x.n;dist=x.dist;return *this;};
};


template <class T=NonSortable>
class NeighList: public std::list<T> {
public:
	NeighList(){};
	virtual ~NeighList(){this->clear();};
	void mysort(){this->sort(compare);};
private:
	static bool compare(T &, T &);
};

template<>
inline bool NeighList<Sortable>::compare(Sortable & x, Sortable & y){
	return x.dist<y.dist;
}
template<>
inline bool NeighList<NonSortable>::compare(NonSortable & x, NonSortable & y){
	return true;
}

typedef NeighList<NonSortable> Neigh;
typedef NeighList<Sortable> Neigh_Vor;

#endif /* NEIGHLISTB_H_ */
