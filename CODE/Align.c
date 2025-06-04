/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * This is part of the DOWSER program
 *
 * DOWSER finds buried water molecules in proteins.
 *  cf. Zhang & Hermans, Proteins: Structure, Function and Genetics 24: 433-438, 1996
 *
 * DOWSER was developed by the Computational Structural Biology Group at the 
 * University of North Carolina, Chapel Hill by Li Zhang, Xinfu Xia, Jan Hermans, 
 * and Dave Cavanaugh.  Revised 1998.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *   filename 
 *      Align: subroutines to rotate a (water) moleucle
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "dowser.h"

#define CALL (void)

extern REAL ScalarProduct3(REAL *a, REAL *b);
extern void CrossProduct(REAL *a, REAL *b, REAL *v);
extern void VectorNormalize3(REAL *a);
void TurnMol(REAL angle, REAL *spindle, REAL *Xyz, int pivot, int Num_atoms);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *     Align (Num_atoms,Xyz,pivot,vector_x,vector)
 *        rotate Xyz about pivot atom,
 *        so that vector_x[] (defined as a function of Xyz) is
 *        parallel with vector[]
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void Align (int Num_atoms, int pivot, REAL *Xyz, REAL *vector_x, REAL *vector)
{
int i;
REAL u[3],v[3],w[3];
REAL angle;

	for (i = 0; i < 3; i++) {
		u[i] = vector_x[i];
		v[i] = vector[i];
	}
	CALL CrossProduct (u, vector, w);

	CALL VectorNormalize3 (u);
	CALL VectorNormalize3 (v);

	angle = acos (ScalarProduct3 (u, v));
	CALL TurnMol (angle, w, Xyz, pivot, Num_atoms);

} /* end Align() */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  void TurnMol (angle, spindle, Xyz, pivot, Num_atoms)
 *       rotate a molecule about a spindle (on atom=pivot, || spindle) by angle
 * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void TurnMol (REAL angle, REAL *spindle, REAL *Xyz, int pivot, int Num_atoms)
{
    REAL u[3], v[3], w[3], s[3];
    REAL cosa1, sina, Pair ();
    REAL ss;
    int i, iat;

    cosa1 = 1. - cos (angle);
    sina = sin (angle);

    /* u is unit spindle parallel to spindle */
	for (i=0;i<3;i++) u[i] = spindle[i];
	CALL VectorNormalize3 (u);

    for (iat = 0; iat < Num_atoms; iat++) {
		for (i=0;i<3;i++) s[i] = *(Xyz+3*iat+i) - *(Xyz+3*pivot+i);

		CALL CrossProduct (u, s, w);
		CALL CrossProduct (u, w, v);

		for (i = 0; i < 3; i++)
			*(Xyz+3*iat+i) += w[i] * sina + v[i] * cosa1;
    }
}               /* end TurnMol() */

#ifdef NOMATH
/* * * * * * * * * * * * * * * * * * * * * * * * *
 *  CrossProduct(a,b,v)
 *     Vector cross product (v = a x b)
 * * * * * * * * * * * * * * * * * * * * * * * * */
void CrossProduct(REAL *a, REAL *b, REAL *v)
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
 REAL ScalarProduct3(REAL *a, REAL *b)
 {
     return(*a* *b + *(a+1)* *(b+1) + *(a+2)* *(b+2));
 }

/* * * * * * * * * * * * * * * * * * * * * * * * *
 * VectorNormalize3(a)
 *     vector a is replaced by a/|a|
 * * * * * * * * * * * * * * * * * * * * * * * * */
void VectorNormalize3(REAL *a)
{
 REAL b;
 b=1./sqrt(*a* *a+ *(a+1)* *(a+1)+ *(a+2)* *(a+2));
 *(a) *= b;
 *(a+1) *= b;
 *(a+2) *= b;
}
#endif
