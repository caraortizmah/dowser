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
 *   Minimize.c
 *     Contains functions for the Fletcher-Reeves-Polak-Ribiere minimization algorithm.
 *     Adapted from Numerical Recipes in C (1990).
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

#define FreeVector free
int ConjGradMinimize(REAL *p, int n, REAL ftol, int *iter, REAL *fret, REAL (*function_ptr)(REAL *, REAL *));
REAL *AllocVector(int nh);
void LineMinimize(REAL p[], REAL xi[], int n, REAL *fret, REAL (*function_ptr)(REAL *, REAL *));
REAL F_atNewVars(REAL x);
void BracketMinimum(REAL *ax, REAL *bx, REAL *cx, REAL *fa, REAL *fb, REAL *fc, REAL (*function_ptr)(REAL));
REAL BrentMinimum(REAL a_input, REAL b_input, REAL c_input, REAL (*function_ptr)(REAL), REAL tol, REAL *x_output);

int NumVars_com=0;
REAL *Vars_com=0,*Grad_com=0;
REAL (*nrfunc)(REAL *, REAL *);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   frprmn
 *     Given a starting point p[0..n-1], Fletcher-Reeves-Polak-Ribiere minimization
 *     is performed on a function func, using its gradient as calculated by a routine
 *     dfunc.  The convergence tolerance on the function value is input as ftol.  
 *     Returned quantities are p (the location of the minimium), iter (the number of
 *     iterations that were performed), and fret (the minimium value of the function).
 *     The routine LineMinimize is called to perform line minimizations.  For more info.
 *     see pp. 317-323.
 *
 *     Returns status code: 0-error 1-okay 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
int ConjGradMinimize(REAL *p, int n, REAL ftol, int *iter, REAL *fret, REAL (*function_ptr)(REAL *, REAL *))
{

#define ITMAX 200
#define EPS 1.0e-10
#define FREEALL FreeVector(xi);FreeVector(h);FreeVector(g);

int j,its;
REAL gg,gam,fp,dgg;
REAL *g,*h,*xi;

g=AllocVector(n);
h=AllocVector(n);
xi=AllocVector(n);
fp=(*function_ptr)(p,xi);
for (j=0;j<n;j++)  {
    g[j] = -xi[j];
	xi[j] = h[j] = g[j];
}
for (its=1;its<=ITMAX;its++)  {
	*iter = its;
	LineMinimize(p,xi,n,fret,function_ptr);
  	/*	if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS))  {  */
	if (fabs(*fret-fp) <= ftol)  {
		FREEALL
		return (1);
	}
	fp = (*function_ptr)(p,xi);
	dgg = gg = 0.0;
	for (j=0;j<n;j++)  {
		gg += g[j]*g[j];
		dgg += (xi[j]+g[j])*xi[j];
	}
	if (gg == 0.0)  {
		FREEALL
		return (1);
	}
	gam = dgg/gg;
	for (j=0;j<n;j++)  {
		g[j] = -xi[j];
		xi[j] = h[j] = g[j]+gam*h[j];
	}
}
#ifdef SHOWERROR
fprintf(stderr,"\n*** Too many iterations in ConjGradMinimize() \n");
#endif
return(0);
} /* end ConjGradMinimize */


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   AllocVector
 *     Allocate a REAL vector of length nh
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

REAL *AllocVector(int nh)
{
REAL *v;
v = (REAL *) malloc ((nh) * sizeof(REAL));
if (!v)  {
	fprintf(stderr,"*** Unable to allocate vector\n");
	exit(1);
}
return v;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  LineMinimize 
 *     Given an n dimensional point p[0..n-1] and an n dimensional direction xi[0..n-1],
 *     moves and resets p to where the function(p) takes on a minimium along the direction
 *     xi from p, and replaces xi by the actual vector displacement that p was moved.
 *     Also returns as fret the value of frun at the returned location p.  This is
 *     actually all accomplished by calling the routines BracketMinimum and BrentMinimum.  For more
 *     info. see pp. 315-317.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void LineMinimize(REAL p[], REAL xi[], int n, REAL *fret, REAL (*function_ptr)(REAL *, REAL *))
{

#ifndef TOL
#define TOL 2.0e-4
#endif


int j;
REAL xx,xmin,fx,fb,fa,bx,ax;
static int allocation=0;

NumVars_com=n;

if (allocation < NumVars_com) {
    FreeVector(Grad_com);
    FreeVector(Vars_com);
    Vars_com=AllocVector(NumVars_com);
    Grad_com=AllocVector(NumVars_com);
    allocation=NumVars_com;
}

nrfunc=function_ptr;
for (j=0;j<n;j++)  {
	Vars_com[j]=p[j];
	Grad_com[j]=xi[j];
}
ax=0.0;
/* xx=1.0; */
xx=1.0e-3;

#define DEBUG
#ifdef DEBUGA
fprintf (stderr,"START MNBRAK routine\n");
#endif
(void) BracketMinimum(&ax,&xx,&bx,&fa,&fx,&fb,F_atNewVars);

#ifdef DEBUGA
fprintf (stderr,"START Brent routine ax= %f, xx= %f, bx= %f, fx= %f, fb= %f\n",ax,xx,bx,fx,fb); 
#endif
*fret=BrentMinimum(ax,xx,bx,F_atNewVars,TOL,&xmin);
#ifdef DEBUGA
fprintf (stderr,"EXIT Brent routine xmin= %f f= %f\n",xmin,*fret);
#endif

for (j=0;j<n;j++)  {
	xi[j] *= xmin;
	p[j] += xi[j];
}
return;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   REAL F_atNewVars(x)
 *     extrapolate p from Vars_com + x * searchvector
 *     and compute function for that value of p
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
REAL F_atNewVars(REAL x)
{

int j;
REAL function_value;
static REAL *Vars_temp;
static int allocation=0;

if (allocation<NumVars_com) {
    FreeVector(Vars_temp);
    Vars_temp = AllocVector(NumVars_com);
}
for (j=0;j<NumVars_com;j++)  Vars_temp[j] = Vars_com[j]+x*Grad_com[j];
function_value = (*nrfunc) (Vars_temp,NULL);
return function_value;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  BracketMinimum
 *    Given a function func, and given distinct initial points ax and bx, this routine
 *    searches in the downhill direction (defined by the function as evaluated at the
 *    initial points) and returns new points ax,bx,cx which bracket a minimum of the
 *    function.  Also returned are the function values at the three points: fa,fb,fc. 
 *    For more info. see pp. 296-298. 
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void BracketMinimum(REAL *ax, REAL *bx, REAL *cx, REAL *fa, REAL *fb, REAL *fc, REAL (*function_ptr)(REAL))
{

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

REAL ulim,u,r,q,fu,dum;

*fa = (*function_ptr)(*ax);
*fb = (*function_ptr)(*bx);
if (*fb > *fa)  {
   SHFT(dum,*ax,*bx,dum); /* this swaps ax and bx */
   SHFT(dum,*fb,*fa,dum); 
}
/* interpolation between a and b */
*cx = (*bx)+GOLD*(*bx-*ax);
*fc = (*function_ptr)(*cx);

while (*fb > *fc)  {
    r = (*bx-*ax)*(*fb-*fc);
    q = (*bx-*cx)*(*fb-*fa);
    u = (*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
    ulim = (*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0)  {
        fu = (*function_ptr)(u);
        if (fu < *fc)  { /* minimum between b and c */
			*ax = (*bx);
			*bx = u;
			*fa = (*fb);
			*fb = fu;
			return;
 		} 
		else if (fu > *fb)  { /* minimum between a and u */
			*cx = u;
			*fc = fu;
			return;
		}
		/* parabolic fit was no good, use default magnification */
		u = (*cx)+GOLD*(*cx-*bx);
		fu = (*function_ptr)(u);
	}
	else if ((*cx-u)*(u-ulim) > 0.0)  { /* parab. fit is between c and its allowed limit */
		fu = (*function_ptr)(u);
		if (fu < *fc)  {
			SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
			SHFT(*fb,*fc,fu,(*function_ptr)(u))
		}
	}
	else if ((u-ulim)*(ulim-*cx) >= 0.0)  { /* lmit parab. u to maximum allowd value */
		u = ulim;
		fu = (*function_ptr)(u);
	}
	else  { /* reject parab. u, use default magnification */
		u = (*cx)+GOLD*(*cx-*bx);
		fu = (*function_ptr)(u);
	}
	SHFT(*ax,*bx,*cx,u) /* eliminate oldest point and continue */
	SHFT(*fa,*fb,*fc,fu)
}

return;
}	


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  BrentMinimum
 *     Given a function f(x), and given three values of x, (ax, bx and cx) for which
 *     ax < bx < cx, and f(bx) < both f(ax) and f(cx), this
 *     routine isolates the minimum to a fractional precision of about tol using
 *     Brent's method.  The location of the minimum is passed back as xmin, and the
 *     function returns the minimum value f(xmin).
 *     For more info. see pp. 299-302.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

REAL BrentMinimum(REAL a_input, REAL b_input, REAL c_input, REAL (*function_ptr)(REAL), REAL tol, REAL *x_output)
{

#define ITMAX2 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10

int iter;
REAL a,b,d,etemp,F_u,F_v,F_w,F_x,p,q,r,tol1,tol2,u,v,w,x,x_mean;
REAL e = 0.0;

a = ((a_input < c_input) ? a_input : c_input);
b = ((a_input > c_input) ? a_input : c_input);
x=w=v=b_input;
F_w=F_v=F_x=(*function_ptr)(x);
for (iter=1;iter<=ITMAX2;iter++)  {
	x_mean = 0.5*(a+b);
	tol2 = 2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-x_mean) <= (tol2-0.5*(b-a)))  {
		*x_output = x;
		return F_x;
	}
	/* construct a parabolic fit */
	if (fabs(e) > tol1)  {
		r = (x-w)*(F_x-F_v);
		q = (x-v)*(F_x-F_w);
		p = (x-v)*q-(x-w)*r;
		q = 2.0*(q-r);
		if (q > 0.0) p = -p;
		q = fabs(q);
		etemp = e;
		e = d;
		/* the following condition determines the acceptability of the parabolic fit;
		   here we tke the golden section step into the larger of the two segments */
		if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))  {
			d = CGOLD*(e=(x >= x_mean ? a-x : b-x));
		}
		else  { /* take the parabolic step */
			d = p/q;
			u = x+d;
			if (u-a < tol2 || b-u < tol2)  d=SIGN(tol1,x_mean-x);
		}
	}
	else {
		d = CGOLD*(e=(x >= x_mean ? a-x : b-x));
    }
	u = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
	F_u = (*function_ptr)(u);
	if (F_u <= F_x)  {
		if (u >= x) a=x;
		else b=x;
		SHFT(v,w,x,u);
		SHFT(F_v,F_w,F_x,F_u);
	}
	else  {
		if (u < x) a=u;
		else b=u;
		if (F_u <= F_w || w == x)  {
			v = w;
			w = u;
			F_v = F_w;
			F_w = F_u;
		}
		else if (F_u <= F_v || v == x || v == w)  {
			v = u;
			F_v = F_u;
		}
	}
}	
	fprintf(stderr,"*** Too many iterations in BRENT\n");
	*x_output = x;
	return F_x;
}
