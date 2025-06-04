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
/* * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  MakeApex() is a subroutine of the qms program
 *    finds the contact points of a fourth sphere with 3 given spheres
 * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <math.h>

#define REAL double

void VectorCopy(REAL *a, REAL *b);
void CrossProduct(REAL *a, REAL *b, REAL *v);
int VectorNormalize3(REAL *a);
REAL ScalarProduct3(REAL *a, REAL *b);
REAL Pair(REAL *a, REAL *b, REAL *u);
REAL DistSq(REAL *a, REAL *b);
int MakeApex(REAL *xyz3fold, REAL *distance, REAL *xyzapex);
int FourBody();

/* * * * * * * * * * * * * * * * * * * * * * * *
 *  MakeApex(xyz3fold,distance,xyzapex)
 *    Find the apex of a tetrahedron given base coordinates
 *        and lengths of the 3 sides
 *    xyz3fold = coordinates of 3 corners of the base
 *    dist = 3 distances where to place the apex
 *    xyzapex = coordinates of 2 apices
 *      needs subroutines from $SIGMA/SIGLIB
 * * * * * * * * * * * * * * * * * * * * * * * */

int MakeApex(REAL *xyz3fold, REAL *distance, REAL *xyzapex)
{
/* input */
REAL *xyz1=xyz3fold,*xyz2=xyz3fold+3,*xyz3=xyz3fold+6;
REAL d1=distance[0],d2=distance[1],d3=distance[2];

/* intermediates */
REAL a,b,c;
REAL xyz4[3],xyz6[3];
REAL cosp, sinp, cosq, sinq;
REAL v12[3],v32[3],v13[3], v54[3],v76[3], v46[3];
REAL e12[3],e32[3],e54[3],e76[3];
REAL xyzn[3];
REAL zvec[3];
REAL d54,d76,u[3];

REAL aa,amax,s;
int j,k,l,l1,k1;

/* output */
REAL *xyzapex1 = xyzapex,*xyzapex2 = xyzapex+3;

/* distances and vectors between the 3 centers */
a = Pair (xyz1,xyz2,v12);
b = Pair (xyz3,xyz2,v32);
c = Pair (xyz1,xyz3,v13);

/* clearly not possible to construct apex */
if (d1+d2 < a || d2+d3 < b || d1+d3 < c) return 1;

/* unit vectors between 3 centers */
VectorCopy (e12,v12); VectorCopy (e32,v32);
if (VectorNormalize3 (e12)) return 1;
if (VectorNormalize3 (e32)) return 1;

/* Construct line number 1 */
/* construct basepoint on 1-2 vector */
cosp = ( d2*d2 + a*a - d1*d1 ) / (2. * a * d2);
for (j=0; j<3; j++) xyz4[j] = xyz2[j] - (d2*cosp) * e12[j];

/* direction of normal from point 5 (in plane) */
CrossProduct (e32,e12,zvec);
if (VectorNormalize3 (zvec)) return 1;
CrossProduct (zvec,e12,e54);
if (VectorNormalize3 (e54)) return 1;

/* Construct line number 2 */
/* construct basepoint on 3-2 vector */
cosp = ( d2*d2 + b*b - d3*d3 ) / (2. * b * d2);
for (j=0; j<3; j++) xyz6[j] = xyz2[j] - (d2*cosp) * e32[j];

/* direction of normal from point 6 (in plane) */
CrossProduct (e32,zvec,e76);
if (VectorNormalize3 (e76)) return 1;

(void) Pair(xyz4,xyz6,v46);
/* construct the projection of the apex on the base */
/* Intersect two lines by solving 2 linear equations */
/* find the best combination of x, y or z components */
amax = 0.;
for (k=0;k<3;k++) {
    l=k+1; if (l==3) l=0;
    aa = e54[k]*e76[l] - e54[l]*e76[k];
    if (aa<0.) aa = -aa;
    if (amax<aa) { amax = aa; k1 = k; l1 = l; }
}
s = ( v46[l1]*e76[k1] - v46[k1]*e76[l1] ) / amax;
for (j=0;j<3;j++) xyzn[j] = xyz4[j] - s*e54[j];

/* check if the distances to the apex exceed the distances to the projection */
if (DistSq (xyz1,xyzn) > d1*d1 || DistSq (xyz2,xyzn) > d2*d2 || DistSq (xyz3,xyzn) > d3*d3)
    return 1;

/* calculate the height of the apex above the base */
s=0.;
for (j=0;j<3;j++) s += (xyz1[j] - xyzn[j])*(xyz1[j] - xyzn[j]);
s = sqrt (d1*d1 - s);

/* construct the apex */
for (j=0;j<3;j++) {
    xyzapex1[j] = xyzn[j] + s * zvec[j];
    xyzapex2[j] = xyzn[j] - s * zvec[j];
}

return 0;
}

/* * * * * * * * * * * * * * * * * * * * * * * *
 * VectorCopy (REAL *a, REAL *b)
 * * * * * * * * * * * * * * * * * * * * * * * */
void VectorCopy (REAL *a, REAL *b)
{ a[0] = b[0]; a[1] = b[1]; a[2] = b[2]; }

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * REAL Pair(a,b,u)
 *      a=coordinates of point 1
 *      b=coordinates of point 2
 *      u=vector b-a
 *      function returns the distance between a and b
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
REAL Pair(REAL *a, REAL *b, REAL *u)
{
int i;
/*
  a,b = 2 positions (input), u=vector a->b (output), pair = |u|
*/
for (i=0;i<3;i++) u[i]=b[i]-a[i];
return (sqrt(ScalarProduct3(u,u)));
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
 * VectorNormalize3(a)
 *     vector a is replaced by a/|a|
 * * * * * * * * * * * * * * * * * * * * * * * * */
#define TOL 1.e-8

int VectorNormalize3(REAL *a)
{
 REAL b;
 b=sqrt(*a* *a+ *(a+1)* *(a+1)+ *(a+2)* *(a+2));
 if (b < TOL) return 1;
 b = 1./b;
 *(a) *= b;
 *(a+1) *= b;
 *(a+2) *= b;
 return 0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * *
*   DistSq: returns distance between 2 points
* * * * * * * * * * * * * * * * * * * * * * * * * */
#undef SQUARE
#define SQUARE(x) (x)*(x)
REAL DistSq(REAL *a, REAL *b)     /* 3-vectors */
 {
 REAL r2;
 int k;

     r2=0.;
     for (k=0;k<3;k++) r2+=SQUARE(*(a+k)-*(b+k));
     return(r2);
}
