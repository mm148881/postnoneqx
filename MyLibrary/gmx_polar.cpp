/*
 * gmx_polar.cpp
 *
 *  Created on: Apr 8, 2012
 *      Author: marchi
 */
#include "gmx_polar.h"
namespace Polar{

int AtomData::nindex=0;
int AtomData::natoms=0;
gmx_bool AtomData::bOXT=false,AtomData::bDipole=false,AtomData::bFilter=false,AtomData::bXplor=false;
gmx_bool AtomData::bTrans=false,AtomData::bFit=false,AtomData::bnoElecPBC=false,
		AtomData::bFluctuations=false;
real * AtomData::w_lst=NULL;
rvec * AtomData::x=NULL;
matrix AtomData::box={1,0,0,0,1,0,0,0,1};
}
