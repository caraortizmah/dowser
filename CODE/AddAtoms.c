/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * This is part of the DOWSER program
 *
 * DOWSER finds buried water molecules in proteins.
 *  cf. Zhang & Hermans, Proteins: Structure, Function and Genetics 24: 433-438, 1996
 *
 * DOWSER was developed by the Computational Structural Biology Group at the 
 * University of North Carolina, Chapel Hill by Li Zhang, Xinfu Xia,Jan Hermans, 
 * and Dave Cavanaugh.  Revised 1998.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 
/* * * * * * * * * * * * * * * * * * * * * * * * * 
 *    S*I*G*M*A - 1996 program
 * revised 1996        J. Hermans, Univ. N. Carolina
 * * * * * * * * * * * * * * * * * * * * * * * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *    Add1Atom(atoms,backback,backat,refatom,atomx,bondlength,bondangle,dihedral)
 *      calculate coordinates of atom atomx as defined by positions of
 *      atoms 1, 2, 3 and standard geometry (dict2)
 *   backback,backat=bond;  ia3=reference;  atomx has no xyz, bonded to backat;
 *     (ia3 may be bonded to backback or to backat)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

#define CALL (void) 
#define SQRT sqrt

REAL ScalarProduct3(REAL *a, REAL *b);
void CrossProduct(REAL *a, REAL *b, REAL *v);
void VectorNormalize3(REAL *a);

int 
    Add1Atom (atoms,backback, backat, refatom, atomx, bondlength, bondangle, dihedral)

int backback, backat, refatom, atomx;
REAL bondlength, bondangle, dihedral;
PDB *atoms;

{
    int i, j, kd, kd1;
    REAL aa, tem, eel;
    REAL a[3], am[9];
    REAL *ex, *ey, *ez;
    REAL tb[3], ta[3];
    int at[4];

    at[0]=backback; at[1]=backat; at[2]=refatom; at[3]=atomx;
    for (i=0;i<4;i++) {
	if (at[i]>=0) continue;
	fprintf (stderr,"REMARK ADD is IMPOSSIBLE.");
	fprintf(stderr," Atoms = %d %d %d %d %d\n",i,backback,backat,refatom,atomx); 
	return 0;
    }
    for (i=0;i<3;i++) {
	if (atoms[at[i]].ftNote) continue;
	fprintf(stderr,"REMARK NO ADD because not all reference xyz set");
	fprintf(stderr," Atoms = %d %d %d %d %d\n",i,backback,backat,refatom,atomx); 
	return 0;
    }
    /*  ****  ex,ey,ez are the local unit vectors */
    /* ex is along bond 1->2; ey is in plane 1,2,3 */
    ex = am;
    ey = am + 3;
    ez = am + 6;

    ex[0] = atoms[backat].XX - atoms[backback].XX;
    ex[1] = atoms[backat].YY - atoms[backback].YY;
    ex[2] = atoms[backat].ZZ - atoms[backback].ZZ;

    /* first ey is along bond 1->3 (except in special case) */
    ta[0] = atoms[refatom].XX;
    tb[0] = atoms[backback].XX;
    ta[1] = atoms[refatom].YY;
    tb[1] = atoms[backback].YY;
    ta[2] = atoms[refatom].ZZ;
    tb[2] = atoms[backback].ZZ;

    for (i=0;i<3;i++) ey[i] = ta[i] - tb[i];

    CALL CrossProduct (ex, ey, ez);	/* ez gets proper direction */

    eel = 0.;
    for (i = 0; i < 3; i++) eel += ez[i] * ez[i];

    if (eel != 0.) {
	CALL VectorNormalize3 (ex);
	CALL VectorNormalize3 (ez);
	CALL CrossProduct (ez, ex, ey);		/* ex and ez normalized and ey recalculated */

	/* the matrix consisting of the three local unit vectors converts */
	/* global xyz to local xyz. Transpose it to reverse its action */
#define TRANSPOSE

#ifdef TRANSPOSE
	tem = am[1];
	am[1] = am[3];
	am[3] = tem;
	tem = am[2];
	am[2] = am[6];
	am[6] = tem;
	tem = am[5];
	am[5] = am[7];
	am[7] = tem;
#endif

	a[0] = -bondlength * cos (bondangle);
	aa = bondlength * sin (bondangle);
	a[1] = aa * cos (dihedral);
	a[2] = aa * sin (dihedral);

	/* transform from local to global coord system */

		atoms[atomx].XX = atoms[backat].XX  + ScalarProduct3 (am     , a);
		atoms[atomx].YY = atoms[backat].YY  + ScalarProduct3 (am + 3 , a);
		atoms[atomx].ZZ = atoms[backat].ZZ  + ScalarProduct3 (am + 6 , a);

    }
    else {
	fprintf (stderr, "REMARK ADDA: zero length vector at atom %d-\n", atomx);
	return 0;
    }

    return 1;

}				/* end Add1Atom() */
/* * * * * * * * * * * * * * * * * * * * * * * * *
 *  CrossProduct(a,b,v)
 *     Vector cross product (v = a x b)
 * * * * * * * * * * * * * * * * * * * * * * * * */
void CrossProduct(a,b,v)
REAL *a,*b,*v;
{
	int i,j,k;
	for (i=0;i<3;i++) {
		j=i+1 ; if (j==3) j=0; k=3-i-j;
		*(v+i)= *(a+j)* *(b+k) - *(a+k)* *(b+j);
	}
}
/* * * * * * * * * * * * * * * * * * * * * * * * *
* REAL ScalarProduct3(a,b)
*     Returns vector scalar product (a . b)
* * * * * * * * * * * * * * * * * * * * * * * * */
REAL ScalarProduct3(a,b)
REAL *a,*b;
{
	return(*a* *b + *(a+1)* *(b+1) + *(a+2)* *(b+2));
}
/* * * * * * * * * * * * * * * * * * * * * * * * *
* VectorNormalize3(a)
*     vector a is replaced by a/|a|
* * * * * * * * * * * * * * * * * * * * * * * * */
void VectorNormalize3(a)
REAL *a;
{
	REAL b;
	b=1./SQRT(*a* *a+ *(a+1)* *(a+1)+ *(a+2)* *(a+2));
	*(a) *= b;
	*(a+1) *= b;
	*(a+2) *= b;
}

