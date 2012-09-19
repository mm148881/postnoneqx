/*
 * Mymap_Spec.hpp
 *
 *  Created on: Oct 29, 2011
 *      Author: second
 */

#ifndef MYMAP_SPEC_HPP_
#define MYMAP_SPEC_HPP_

int Tier::mx=0;
int Tier::my=0;
int Tier::mz=0;


template <>
void Mymap<int,int>::insert(std::pair<int,int> b){
	  pair<std::map<int,int>::iterator,bool> ret;
	  ret=std::map<int,int>::insert(b);
	  if (ret.second==false) ret.first->second+=b.second;
}
template <>
void Mymap<Tier,int>::insert(std::pair<Tier,int> b){
	  pair<std::map<Tier,int>::iterator,bool> ret;
	  ret=std::map<Tier,int>::insert(b);
	  if (ret.second==false) ret.first->second+=b.second;
}


#endif /* MYMAP_SPEC_HPP_ */
