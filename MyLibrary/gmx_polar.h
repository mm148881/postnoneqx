/*
 * gmx_charges.h
 *
 *  Created on: May 25, 2011
 *      Author: marchi
 */

#ifndef GMX_CHARGES_H_
#define GMX_CHARGES_H_

#include "copyrite.h"
#include "filenm.h"
#include "macros.h"
#include "pbc.h"
#include "smalloc.h"
#include "statutil.h"
#include "vec.h"
#include "xvgr.h"
#include "typedefs.h"

#include "nbsearch.h"
#include "trajana.h"
#include "futil.h"
#include "index.h"
#include "strdb.h"
#include "string2.h"
#include "pbc.h"
#include "rmpbc.h"
//#include "mkl.h"
#include "gmx_fatal.h"
#include "do_fit.h"
#include "types/block.h"
#include "AtomIndex.h"
#include "Charges.h"
#include "AtomsPDB.h"
#include "Phi.h"
#include "Dielectric.h"
#include "Dipoles.h"
#include "RhoDip.h"
#include "Averages.h"
namespace Polar {
struct AtomData{
	static int nindex;
	static int natoms;
	static gmx_bool bOXT,bDipole,bFilter,bXplor,bTrans,bFit,bnoElecPBC,bFluctuations;
	static real * w_lst;
	static matrix box;
	static rvec * x;
};
static matrix MyCO={{0,0,0},{0,0,0},{0,0,0}};
static rvec supercell={6,6,6};

static AtomIndex * cidx;
static AtomIndex cidx_ref;
static atom_id * cindex;

enum { euSel,euRect, euTric, euCompact, euNR};
enum ElecType {Slt, Sol, Tot};
static rvec * xrefs;
static int NatomRef=0;
static rvec ngrid={0,0,0};
static rvec  trans = {0,0,0};
static Charges * chg;
static AtomsPDB * atPDB=NULL;

static int  ngrps=0;
static Dipoles MyDip0;
static Dipoles MyDip;
static t_block MyMols;
static void center_x(int ecenter,rvec x[],matrix box,int n,int nc,atom_id ci[]){
    int i,m,ai;
    rvec cmin,cmax,box_center,dx;

    if (nc > 0) {
        copy_rvec(x[ci[0]],cmin);
        copy_rvec(x[ci[0]],cmax);
        for(i=0; i<nc; i++) {
            ai=ci[i];
            for(m=0; m<DIM; m++) {
                if (x[ai][m] < cmin[m])
                    cmin[m]=x[ai][m];
                else if (x[ai][m] > cmax[m])
                    cmax[m]=x[ai][m];
            }
        }
        calc_box_center(ecenter,box,box_center);
        for(m=0; m<DIM; m++)
            dx[m] = box_center[m]-(cmin[m]+cmax[m])*0.5;

        for(i=0; i<n; i++)
            rvec_inc(x[i],dx);
    }
}

static void put_molecule_com_in_box(int unitcell_enum,int ecenter,
                                    t_block *mols,
                                    int natoms,t_atom atom[],
                                    int ePBC,matrix box,rvec x[])
{
    atom_id i,j;
    int     d;
    rvec    com,new_com,shift,box_center;
    real    m;
    double  mtot;
    t_pbc   pbc;
    calc_box_center(ecenter,box,box_center);
    set_pbc(&pbc,ePBC,box);
    if (mols->nr <= 0)
        gmx_fatal(FARGS,"There are no molecule descriptions. I need a tpr file for this pbc option.");
    for(i=0; (i<mols->nr); i++) {
        /* calc COM */
        clear_rvec(com);
        mtot = 0;
        for(j=mols->index[i]; (j<mols->index[i+1] && j<natoms); j++) {
            m = atom[j].m;
            for(d=0; d<DIM; d++)
                com[d] += m*x[j][d];
            mtot += m;
        }
        /* calculate final COM */
        svmul(1.0/mtot, com, com);
	for(d=0; d<DIM; d++)
        /* check if COM is outside box */
        copy_rvec(com,new_com);
        switch (unitcell_enum) {
        case euRect:
            put_atoms_in_box(box,1,&new_com);
            break;
        case euTric:
            put_atoms_in_triclinic_unitcell(ecenter,box,1,&new_com);
            break;
        case euCompact:
            put_atoms_in_compact_unitcell(ePBC,ecenter,box,1,&new_com);
            break;
        }
        if (debug) fprintf(debug," My debug %d %d \n",ePBC,ecenter);
        rvec_sub(new_com,com,shift);
        if (norm2(shift) > 0) {
            if (debug)
                fprintf (debug,"\nShifting position of molecule %d "
                         "by %8.3f  %8.3f  %8.3f\n", i+1, PR_VEC(shift));
            for(j=mols->index[i]; (j<mols->index[i+1] && j<natoms); j++) {
                rvec_inc(x[j],shift);
            }
        }else if(debug) fprintf (debug,"\n My Debug 2 %d  %f %f %f \n", i+1, new_com[XX], new_com[YY], new_com[ZZ]);
    }
}

/*! \brief
 * Template analysis data structure.
 */

}
#endif /* GMX_CHARGES_H_ */
