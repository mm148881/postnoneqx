/*
 * gmx_charges.h
 *
 *  Created on: May 25, 2011
 *      Author: marchi
 */

#ifndef GMX_VORONOI_H_
#define GMX_VORONOI_H_


#include "AtomIndex.h"
#include "Atoms.h"
#include <iostream>
#include <fstream>
#include "Rho.h"

namespace gmx_Voronoi{
static Atoms atm;
static std::ofstream VoroOut;
static int iGrid[DIM]={-1,-1,-1};
static RhoAverage * MyRho0;
static Charges * chg;
static AtomIndex * cidx;
static AtomIndex cidx_ref;
static int nindex;
static atom_id * cindex;

static int Tot;
enum ElecType {Slt=0};
static rvec  trans = {0,0,0};
static int  ngrps=0;
static matrix MyCO={{0,0,0},{0,0,0},{0,0,0}};
static int natoms=0;
static rvec ngrid={0,0,0};
static AtomsPDB * atPDB=NULL;
}
#endif
