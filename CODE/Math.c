/* * * * * * * * * * * * * * * * * * * * * * * * * 
 *    S*I*G*M*A - 1996 program
 * revised 1996        J. Hermans, Univ. N. Carolina
 *    adapted to DOWSER
 * * * * * * * * * * * * * * * * * * * * * * * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * 
 *   mathematical routines:
 *
 *    ScalarProduct3	=  scalar product
 *    CrossProduct	=  vector product
 *    VectorNormalize3	=  normalize a vector
 *    MatrixInv3	=  invert 3x3 matrix
 *    MatrixXVector3	=  matrix times 3-vector
 *    MatrixInv		=  matrix inverter
 *    LeastCommonMultiple
 *    DistSq		= distance**2 
 *    Pair(a,b,u)
 * Tamar Schlick's random number routines:
 * REAL RandU () 	= uniform distribution, single value
 * RandU_vec (n, vec) 	= uniform distribution, array of values
 * int GetRandomSeed () = returns value of current "seed"
 * PutRandomSeed ()     = changes value of curent "seed"
 * RandN_vec1 (n, vec)  = normal distribution, array of values
 * RandN_vec2 (n, vec, mean, var)
 * RandN_vec3 (n, vec, mean, var)
 * RandN_3 ()           = normal distribution, 3 values.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  malloc(): Storage to invert matrix of arbitrary size in MatrixInv().
 *  malloc(): free() and again malloc() if need more memory than before.
 * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "dowser.h"

int Random_seed; /* external variable */

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

/* * * * * * * * * * * * * * * * * * * * * * * * * 
 * int MatrixInv3(r,N,deter)
 *     vector a is replaced by a/|a|
 *    r is the 3 x 3 matrix that is replaced by its inverse
 *    the function returns FALSE if OK, TRUE if singular matrix
 * * * * * * * * * * * * * * * * * * * * * * * * */
int MatrixInv3(REAL *r, int N, REAL *deter)
{
 REAL cofa1,cofa2,cofa3,cofa4,cofa5,cofa6,cofa7,cofa8,cofa9;
 REAL rdet; /* reciprocal of determinant */
/*     compute cofactors */
#define RR(j,i) (*(r+3*i+j))
 cofa1=RR(1,1)*RR(2,2)-RR(1,2)*RR(2,1);
 cofa2=RR(1,2)*RR(2,0)-RR(1,0)*RR(2,2);
 cofa3=RR(1,0)*RR(2,1)-RR(1,1)*RR(2,0);
 cofa4=RR(0,0)*RR(2,2)-RR(0,2)*RR(2,0);
 cofa5=RR(0,1)*RR(2,0)-RR(0,0)*RR(2,1);
 cofa6=RR(0,0)*RR(1,1)-RR(0,1)*RR(1,0);
 cofa7=RR(2,1)*RR(0,2)-RR(0,1)*RR(2,2);
 cofa8=RR(0,1)*RR(1,2)-RR(1,1)*RR(0,2);
 cofa9=RR(1,0)*RR(0,2)-RR(0,0)*RR(1,2);
/*     compute determinant  */
 *deter=RR(0,0)*cofa1+RR(0,1)*cofa2+RR(0,2)*cofa3;
 if (*deter<.0001) {
     fprintf(stderr,"INV3x3:******singular matrix; not inverted*****\n");
     return(TRUE);
 }
 /*   compute inverse matrix */
 rdet=1./(*deter);
 RR(0,0)=cofa1*rdet ; RR(1,0)=cofa2*rdet ; RR(2,0)=cofa3*rdet; 
 RR(1,1)=cofa4*rdet ; RR(2,1)=cofa5*rdet ; RR(2,2)=cofa6*rdet;
 RR(0,1)=cofa7*rdet ;RR(0,2)=cofa8*rdet ; RR(1,2)=cofa9*rdet ;
 return(FALSE);

} /* end inv3x3 */

/* * * * * * * * * * * * * * * * * * * * * * * * * *
 * MatrixXVector3(a,b,c)
 *    3x3 matrix (a) times 3-vector (b)
 *    a x b => c 
 * * * * * * * * * * * * * * * * * * * * * * * * * */

void MatrixXVector3(REAL *a, REAL *b, REAL *c)
{
 int i,j;
 for (i=0;i<3;i++) {
     *(c+i)=0.0 ;
     for (j=0;j<3;j++) *(c+i)+= *(a+j+3*i)* *(b+j);
 }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * *
 *   dist_sq: returns distance between 2 points
 * * * * * * * * * * * * * * * * * * * * * * * * * */
REAL DistSq(REAL *a, REAL *b)	/* 3-vectors */
{
    REAL r2;
    int k;
    REAL aaa; 

    r2=0.;
    for (k=0;k<3;k++) {
	aaa = (*(a+k)-*(b+k));
	r2+=aaa*aaa;
    }
    return(r2);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * REAL Pair(a,b,u)
 *	a=coordinates of point 1
 *	b=coordinates of point 2
 *	u=vector b-a
 *	function returns the distance between a and b
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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * MatrixInv(A,N,d)
 *
 *	minv inverts the N*N matrix A, using the standard Gauss-
 *	Jordan method. The determinant is also calculated.
 *
 *	a[N*N) = matrix to be inverted; delivered with the inverse
 *	N = order of matrix a
 *	d = delivered with the determinant of A
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define TRUE 1
#define FLOAT double
int MatrixInv(REAL *A, int N, REAL *determinant)
{ 
    static int old_N=0;
    static FLOAT *coef,*consta;
    FLOAT *acoef,*aconst,*const_sub,*coef_sub,*top_const,*top_coef;
    FLOAT *coef_from,*coef_to,*const_from,*const_to;
    REAL *anA,*lll;
    int N_sub,row,col,row_sub,col_sub;
    FLOAT maxi_abs_coef,maxi_coef,aaa,lead_const,lead_coef,s;
    int i,row_maxi,col_maxi,row_scan,col_scan;


    if (old_N!=N) {
	if (old_N!=0) free(coef);
	coef=(FLOAT *) malloc(N*N*2*sizeof(FLOAT));
	consta=coef+N*N;
    }

    /* initialize the 2 allocated arrays */
    anA=A; acoef=coef; aconst=consta;
    for (row=0;row<N;row++) {
    for (col=0;col<N;col++) {
	*(acoef++)= *(anA++);
	*(aconst++)=(row==col)?1.:0.; 
    }}

    /* REDUCTION SECTION */
    *determinant=1.;
    for (row=0;row<N;row++) {
	col=row;

	/* find largest coefficient in sub-array */
	coef_sub=coef+row;
	maxi_coef= *(coef_sub+N*row);
	maxi_abs_coef=(maxi_coef>0.)?maxi_coef:-maxi_coef;
	row_maxi=row; col_maxi=row;
	for (row_scan=row;row_scan<N;row_scan++) {
	    aaa= *(acoef=coef_sub+N*row_scan);
	    if (aaa<0.) aaa= -aaa;
	    for (col_scan=col;col_scan<N;col_scan++) {
		if (aaa>maxi_abs_coef) {
		    maxi_abs_coef= aaa;
		    maxi_coef= *acoef;
		    row_maxi=row_scan; col_maxi=col_scan;
		}
		acoef++;
	    }
	}
	(*determinant)*=maxi_coef;
	if ((*determinant)==0.) return 1; /* singularity */

	/* SWAP of rows (row) and (row_maxi) */
	coef_from=coef+row*N; coef_to=coef+row_maxi*N;
	const_from=consta+row*N; const_to=consta+row_maxi*N;
	for (i=0;i<N;i++) {
	    aaa= *(coef_from); *(coef_from++)= *coef_to; *(coef_to++)=aaa;
	    aaa= *(const_from); *(const_from++)= *const_to; *(const_to++)=aaa;
	}
	/* SWAP of columns (col=row) and (col_maxi) */
	coef_from=coef+col; coef_to=coef+col_maxi;
	const_from=consta+col; const_to=consta+col_maxi;
	for (i=0;i<N;i++) {
	    aaa= *(coef_from); *(coef_from)= *coef_to; *(coef_to)=aaa;
	    aaa= *(const_from); *(const_from)= *const_to; *(const_to)=aaa;
	    coef_from+=N; coef_to+=N; const_from+=N; const_to+=N;
	}

	const_sub=consta+row*(N+1);
	coef_sub=coef+row*(N+1);
	N_sub=N-row;

	/* divide entire const row (row) by the lead coef (maxi_coef)
	   and subtract from each lower element the product of
	   this corrected const * (col=row,row=same lower row) coeff */
	top_const=consta+row*N;
	for (col_sub=0;col_sub<N;col_sub++) {
	    lead_const= ((*top_const)/=maxi_coef);
	    for (row_sub=1;row_sub<N_sub;row_sub++) {
		(*(top_const+N*row_sub)) -= (*(coef_sub+N*row_sub))*lead_const;
	    }
	    top_const++;
	}

	/* divide partial row (row) by the leading coefficient  (maxi_coef) 
	   and subtract from all lower coefs in the product of this
	   corrected coeff* leading coeff in that row */
	top_coef=coef_sub+1; 
	for (col_sub=row+1;col_sub<N;col_sub++) {
	    lead_coef=((*top_coef)/=maxi_coef);

	    /* subtract from the following rows
		product of (head of subrow)(head of subcol) */
	    for (row_sub=1;row_sub<N_sub;row_sub++) {
		(*(top_coef+N*row_sub)) -= (*(coef_sub+N*row_sub))*lead_coef;
	    }
	    top_coef++; 
	}
    } /* end REDUCTION */

    /* BACK SUBSTITUTION */
    for (row=N-1;row>=0;row--) {
	coef_sub=coef+row*N;
	const_sub=consta+row*N;
	for (col=0;col<N;col++) { /* in all the constant columns */
	    s = (*(aconst=const_sub+col));
	    for (col_sub=row+1;col_sub<N;col_sub++) { /* sum c*x as far known */
		aconst+=N; /* step down coef col */
		acoef=coef_sub+col_sub; /* step along coef row */
		s -= (*aconst) * (*acoef);
	    }
	    *(const_sub+col)=s;
	}
    }

#ifdef DEBUG
    fprintf(stderr,"\n");
    for (row=0;row<N;row++) {
	for (col=0;col<N;col++) {
	    fprintf(stderr,"%10.3f ",consta[col+N*row]);
	}
	fprintf(stderr,"\n");
    }
#endif

    /* place the result in the array A */
    for (i=0;i<N*N;i++) *(A+i)=(REAL) *(consta+i);
    return 0;
} /* end minv() */

/* * * * * * * * * * * * * * * * * * * *
 *  int LeastCommonMultiple(j,k)
 *    returns the least common multiple of 2 integers.
 *    the range of prime factors is limited by the table.
 *    if table is too small, the function returns 0.
 *    if one of the args is zero, then also return 0.
 *    if an arg is negative, get positive answer.
 * * * * * * * * * * * * * * * * * * * */
int LeastCommonMultiple(int j, int k)
{
#define NUM_PRIMES 10
    static int primes[NUM_PRIMES]={2,3,5,7,11,13,17,19,25,29};
    int i,maxi,mini,ix,inx,ct,maxa;

    if (j*k==0) return (0);
    if (j>k) { maxi=j; mini=k; }
    else { maxi=k; mini=j; }
    maxa=maxi;
    /* find prime factors in the smallest */
    for (inx=0;inx<NUM_PRIMES;inx++) {
	ix=primes[inx]; ct=0;
	while (!(mini%ix)) { mini/=ix; ct++;}

	/* divide the largest by the prime factors */
	while (ct>0 && !(maxa%ix)) { maxa/=ix; ct--; }

	/* multiply the largest by the remaining ones */
	for (i=0;i<ct;i++) maxi*=ix; /* remaining ones */
    }
    if (mini>1) return(0);
    else return (maxi);
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifdef OLDRANDU
REAL RandU (ix)
long int *ix;
{
long int iy;
REAL yfl;
static long int masker=2147483647;

    iy= *ix * 452807053;
    iy=iy & masker;
    yfl=(REAL) iy;
    yfl *= .4656612875e-9;
    *ix = iy+1;
    return( yfl );
}
#else

#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  Copyright (c) 1991 by Tamar Schlick and Samuel A. Figueroa
 *  (c-version jhermans, Nov 94)
 * Following are subprograms RandU, RandU_vec, GetRandomSeed, and PutRandomSeed
 * - all related to the generation of uniformly distributed random numbers -
 * and routines RandN_3, RandN_vec1, RandN_vec2 - related to generation of
 * normally-distributed random numbers.  The method used for the 
 * uniform distribution is based on CACM vol. 31 no. 10 (Oct. 1988),
 * pp. 1192-1201.  This article asserts that no random
 * number generator other than the one implemented here should be used
 * unless it can be proven to be better than this one, which throughout
 * the years has proven to be quite good.  These routines work correctly
 * on VAX's running under the VMS operating system.  THEY ARE NOT GUARAN-
 * TEED TO WORK ON ANY OTHER COMPUTER.  The article mentioned above des-
 * cribes how to implement this method on a wide variety of computers.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define A 16807
#define m 2147483647
#define q 127773
#define r 2836
#ifndef REAL
#define REAL float
#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The function RandU returns a pseudo-random number from a uniform dis-
 * tribution between 0 and 1 exclusive.  The subroutine PutRandomSeed (see be-
 * low) should be called (once) before the first time this function is
 * called.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
REAL RandU (void)
{
    static REAL minv = ((REAL) 1) / m;

    Random_seed = A * (Random_seed % q) - r * (Random_seed/q);
    if (Random_seed <= 0) Random_seed += m;
    return  Random_seed * minv;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The subroutine RandU_vec generates a vector of n pseudo-random numbers,
 * each of which is from a uniform distribution between 0 and 1 exclu-
 * sive, storing this vector in the parameter vec.  The subroutine
 * PutRandomSeed (see below) should be called (once) before the first time
 * this subroutine is called.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void RandU_vec (int n, REAL *vec) 
{
    static REAL minv = ((REAL) 1) / m;
    REAL *avec;
    int i;

    if (n < 1) return;
    avec=vec;
    for (i=0;i<n;i++) {
	Random_seed = A * (Random_seed % q) - r * (Random_seed/q);
	if (Random_seed <= 0) Random_seed += m;
	*(avec++) = Random_seed * minv;
    }
    return;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The subroutine GetRandomSeed returns the current
 * value of the Random_seed used in RandU and RandU_vec
 * to generate random numbers.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int GetRandomSeed (void)
{
    return Random_seed;
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The subroutine PutRandomSeed sets the Random_seed (used in RandU and RandU_vec
 * to generate random numbers) to the value of its (only) parameter iseed, pro-
 * vided iseed is between 0 and the modulus m (see the parameter state-
 * ment below) exclusive.  Otherwise, PutRandomSeed sets the Random_seed to 1.  This
 * subroutine should be called (once) before the first time RandU or RandU_vec
 * is called.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void PutRandomSeed (int iseed)
{
      if (iseed < 0 && iseed < m) Random_seed = iseed;
      else Random_seed = 1;
      return;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The subroutine RandN_vec1 generates a vector of n pseudo-random numbers,
 * each of which is from a standard normal distribution (mean zero, vari-
 * ance one), based on the method of Odeh and Evans: "The Percentage
 * Points of the Normal Distribution," Applied Statistics, 23 (1974),
 * pp. 96-97.  For a mean mu other than zero and/or a variance vector
 * sigma other than one, set:     vec(i) = mu + sqrt(sigma(i))*vec(i).
 * See also the subroutine RandN_vec2 below.  The subroutine PutRandomSeed (above)
 * should be called once before the first time this subroutine is called.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void RandN_vec1 (int n, REAL *vec)
#define p0 (REAL) -.322232431088e0
#define p1 (REAL) -1e0
#define p2 (REAL) -.342242088547e0
#define p3 (REAL) -.204231210245e-1
#define p4 (REAL) -.453642210148e-4
#define q0 (REAL) .99348462606e-1
#define q1 (REAL) .588581570495e0
#define q2 (REAL) .531103462366e0
#define q3 (REAL) .10353775285e0
#define q4 (REAL) .38560700634e-2

{
    REAL temp,*avec;
    int i;
    if (n < 1) return;
    RandU_vec(n, vec); avec=vec;
    for (i=0;i<n;i++) {
	temp = *avec;
	if (temp > 0.5) *avec = (REAL) 1 - *avec;
	*avec = sqrt(log(((REAL) 1)/(*avec * *avec)));
	*avec +=
	    ((((*avec * p4 + p3) * *avec + p2) *
	    *avec + p1) * *avec + p0) /
	    ((((*avec * q4 + q3) * *avec + q2) *
	    *avec + q1) * *avec + q0);
	if (temp < 0.5) *avec = -*avec;
	avec++;
    }
    return;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The subroutine RandN_vec2 generates a vector of n pseudo-random numbers,
 * the ith component of which is from a normal distribution with the
 * specified mean and variance var[i].  This vector is stored in the
 * parameter vec.  The subroutine PutRandomSeed (see above) should be called
 * (once) before the first time this subroutine is called.  This subrou-
 * tine offers an alternative to RandN_vec1 above.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void RandN_vec2 (int n, REAL *vec, REAL mean, REAL *var)
{
      int i;
      REAL *avec, *avar;
      
      RandN_vec1(n, vec);
      avec=vec;
      avar=var;
      for (i=0;i<n;i++) *(avec++) *= sqrt (*avar);

      if (mean != 0) { avec=vec; for (i=0;i<n;i++) *(avec++) += mean; }
      return;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The routine RandN_vec3 assigns to each component i of the vector vec[n] a
 * random number from a Gaussian distribution with prescribed mean and
 * variance var[i]. vec[i] is assigned as a weighted average of m irv's
 * {r1, r2, r3, ... rm} chosen from a uniform dist. between -1 and 1:
 * vec[i]  <--  alpham * (r1 + r2 + .... + rm), where alpham is chosen
 * so that vec[i] (which has a Gaussian dist. by the central limit thm.
 * with mean 0) will have a variance satisfying  < vec[i]**2 > = var[i].
 * From   < vec[i]**2 >  = alpham **2 *  < (r1 + r2 + ... + rm) **2 >
 *                       = alpham **2 * m * 1/3         =  var[i],
 * alpham  = sqrt ( 3 * var[i] / m  ). Since the sum of irv's with
 * normal distributions has a normal density with a mean and variance
 * that are the sums of the individual means and variances, 
 *  < sum {vec[i]**2} >  = sum {var[i]}.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void RandN_vec3 (int n, REAL mean, REAL *vec, REAL *var)
#define zero (REAL) 0.0
#define one (REAL) 1.0
#define two (REAL) 2.0
#define three (REAL) 3.0
#define mlge (REAL) 20
{
    int i,j;
    REAL sumran,alpham,*avar,*avec;
    avar=var; avec=vec;
    for (i=0;i<n;i++) {
         sumran = zero;
	 for (j=0;j<mlge;j++)  sumran += one - two *  RandU();
         alpham = sqrt (  (*(avar++) * three) / (REAL) (mlge)  );
         *avec = alpham * sumran;
    }
    if (mean != zero) {
	avec=vec; for (i=0;i<n;i++) *(avec++) += mean;
    }
    return;
}

#define NEG_TWO (REAL) -2.0
#define POS_TWO (REAL) 2.0

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The function RandN_3 returns a vector of 3 pseudo-random numbers from 
 * a standard normal distribution (mean zero, variance one)
 * A saving is obtained by using one "left-over" item the next time the
 * function is called.
 *   This routine uses a different method from the above three.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void RandN_3 (REAL* ran_force) 

{
static REAL a, b;
static int flipflop = 0;
      
    if (flipflop == 0) {
	a = sqrt(NEG_TWO * log(RandU()));
	b = POS_TWO * PI * RandU();
	flipflop = 1;
	ran_force[0] = a * cos(b);
	ran_force[1] = a * sin(b);
	a = sqrt(NEG_TWO * log(RandU()));
	b = POS_TWO * PI * RandU();
	ran_force[2] = a * cos(b);
    } else {
	flipflop = 0;
	ran_force[0] = a * sin(b);
	a = sqrt(NEG_TWO * log(RandU()));
	b = POS_TWO * PI * RandU();
	ran_force[1] = a * cos(b);
	ran_force[2] = a * sin(b);
    }
}
