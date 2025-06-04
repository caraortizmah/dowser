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
 *  qms program: find contact points of molecular surface
 *    This replaces the Connolly ms program, which requires much more storage
 *        but ONLY to find the internal cavities.
 *    Usage:
 *    qms ATOMFILE #atoms RADIUSFILE #radii
 * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <math.h>
#include "dowser.h"

#undef REAL
#define REAL double
#define SHORT short int
#define TRUE 1
#define FALSE 0
#define BOXSIZE 2
#define NUMPAIRMAXINIT 200
#define MAXCHAR 256
#define MAXRADII 100

#include "pair.H"

REAL *Xyz;
int Num_atoms;
void HashXyz_l(void);
void ScanHash(REAL *xcop, int *pairs, int *npair, int Num_box_cutoff);
extern void initVDW(char *filename, int num_atoms, PDB *atoms);
extern REAL DistSq(REAL r1[3], REAL r2[3]);
extern int MakeApex(REAL *xyz3fold, REAL *distance, REAL *xyzapex);

int main(int argc, char *argv[])
{

#undef DEBUG

int *atomtype;
REAL *aXyz, *radius, *ms_rad; 

REAL aa;
REAL xyztriangle[9], distance[3], xyzapex[6];

int i,j,atom0,atom1,atom2,atom3,i1,i2,i3,apex;
int numpair;
static int numpairmax=NUMPAIRMAXINIT;
int *pairs, *apair;
int surfcount=0;
int Num_radii;
FILE *inputfile; char record[MAXCHAR];
float x1,x2,x3;
int myerror = FALSE;
REAL largest_radius; int number_of_boxes_scan;
REAL probe_radius;

if (argc<6) {
    fprintf(stderr,"\nCORRECT USAGE is: qms ATOMFILE #atoms RADIUSFILE #radii probe_radius\n\n");
    exit (1);
}
sscanf (argv[2],"%d",&Num_atoms);
    Xyz = (REAL *) malloc (Num_atoms*3*sizeof(REAL));
    radius = (REAL *) malloc (Num_atoms*sizeof(REAL));
sscanf (argv[4],"%d",&Num_radii);
    if (Num_radii<MAXRADII) Num_radii=MAXRADII;
    ms_rad = (REAL *) malloc (Num_radii*sizeof(REAL));
sscanf (argv[5],"%f",&x1);
    probe_radius = x1;

/*  ms.dow format:
  31.75800  13.35800 -13.67300    1    2    0 ARG N
    ms.rad format
    1   2.40000
*/

    /* read ms radii from ms.rad */
    if (!(inputfile=fopen(argv[3],"r"))) {
	fprintf (stderr,"Cannot open %s\n",argv[1]);
	exit (1);
    }
    while (TRUE) {
	fgets(record, MAXCHAR, inputfile);
	if( feof(inputfile) ) break;
	sscanf(record, "%d %f", &j, &x1);
	if (j>Num_radii) {
	    fprintf (stderr,"Not enough space allocated (%d) for ms radii, need at least %d\n",Num_radii,j);
	    myerror = TRUE;
	}
	else ms_rad[j-1] = x1;
    }
    if (myerror) exit (1);

    /* read atom data from ms.dow */
    if (!(inputfile=fopen(argv[1],"r"))) {
	fprintf (stderr,"Cannot open %s\n",argv[1]);
	exit (1);
    }
    atom0=0;
    aXyz = Xyz; 
    while (TRUE) {
	fgets(record, MAXCHAR, inputfile);
	if( feof(inputfile) ) break;
	sscanf(record, "%f %f %f %d", &x1, &x2, &x3, &j);
	*(aXyz++) = x1;
	*(aXyz++) = x2;
	*(aXyz++) = x3;
	radius[atom0++] = ms_rad[j-1] + probe_radius;
    }

    largest_radius = 0.;
    for (atom0=0; atom0<Num_atoms; atom0++)
	if (largest_radius<radius[atom0]) largest_radius = radius[atom0];
      number_of_boxes_scan = (int) (2.*largest_radius/BOXSIZE + .5);

    HashXyz_l();
    /* *************************** */
    pairs = (int *) malloc (numpairmax * sizeof (int));
    for (atom0=0; atom0<Num_atoms; atom0++ ) {
	aXyz=Xyz+3*atom0;
/*
	fprintf (stderr," coordinates: %6.2f, %6.2f, %6.2f\n",aXyz[0],aXyz[1],aXyz[2]);
*/
	/* find the near neighbors of this atom */
	numpair=0;
	ScanHash (Xyz+3*atom0, pairs, &numpair, number_of_boxes_scan);
/*
	fprintf (stderr," atom %d, %d neighbors\n", atom0, numpair);
	for (j=0;j<numpair;j++) fprintf(stderr," %d", pairs[j]);
	fprintf(stderr,"\n");
*/

	/* transfer the coordinates and the radii */
	for (j=0;j<3;j++) xyztriangle[j] = *(Xyz + 3*atom0 +j);
	distance[0] = radius[atom0];

	/* find second atom from the near neighbors */
	for (i1=0;i1<numpair;i1++) {
	    atom1 = pairs[i1]-1;
	    if (atom1<=atom0) continue;
	    for (j=0;j<3;j++) xyztriangle[3+j] = *(Xyz + 3*atom1 +j);
	    distance[1] = radius[atom1];

	    /* find third atom from the near neighbors */
	    for (i2=0;i2<numpair;i2++) {
		atom2 = pairs[i2]-1;
		if (atom2<=atom0 || atom2<=atom1) continue;
		for (j=0;j<3;j++) xyztriangle[6+j] = *(Xyz + 3*atom2 +j);
		distance[2] = radius[atom2];

		/* locate the two apices */
		if (MakeApex (xyztriangle,distance,xyzapex)) continue;
		else for (apex=0;apex<6;apex+=3) {
		    aXyz = xyzapex+apex;
		    for (i3=0;i3<numpair;i3++) {
			atom3 = pairs[i3]-1;
			if (atom3==atom0 || atom3 == atom1 || atom3==atom2) continue;
			aa = radius[atom3];
			if (DistSq(aXyz,Xyz+3*atom3) < aa*aa) break;  
		    }
		    if (i3==numpair) {
			surfcount++;
#ifdef DEBUG
			fprintf(stderr,"Surface contact atoms %d %d %d, %6.2f %6.2f %6.2f\n",
			    atom0, atom1, atom2, aXyz[0],aXyz[1],aXyz[2]);
#endif
			fprintf(stdout,"    1    0    0 1 %8.3f %8.3f %8.3f 0\n",
			    aXyz[0],aXyz[1],aXyz[2]);
/* format of output xms.dow
,,,,1,,,,0,,,,0,1,,,33.035,,,15.370,,-13.962,0
*/
		    }
		}
	    }
	}
    }

    fprintf(stderr,"Number of surface points = %d\n",surfcount);
    return 0;

} /* end qms() */

#undef DEBUG

#undef XYZ
#define XYZ(i,j) (*(Xyz+3*j+i-4))

REAL xzero,yzero,zzero; /* Set in HashXyz() */

/* * * * * * * * * * * * * * * * * * * * * * * * * 
 *    S*I*G*M*A - 1996 program
 * revised 1996        J. Hermans, Univ. N. Carolina
 * * * * * * * * * * * * * * * * * * * * * * * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * 
 *  HashXyz_l()
 *
 *  create the hashing representation used in naybor (every call)
 *
 *  hashing representation is as follows:
 *  for each atom zbox and nextatom exist;
 *  for each line exist mini-maxes of zbox and
 *  a pointer to the atom in that line with lowest z;
 *  nextatom points to next highest y,z and ==0 for last atom in line.
 *   malloc: Links between atoms and low and high of each line.
 *   malloc: Number lines depends on volume.
 *   malloc: Free() and again malloc() if need more.
 * * * * * * * * * * * * * * * * * * * * * * * * * */
/* #define DEBUG */


void HashXyz_l(void)
{
REAL abox0[3],abox1[3];

int where,lastwhere;
static int new=TRUE;
static int num_lines_old=0;
int bigger,equal;
REAL a2,a3,aa;
int iat,i,ix,iy,iz;
REAL *x_pt;
REAL aboxsz = BOXSIZE;

/* ***  determine number of boxes from the atomic coordinates */
for (i=1;i<=3;i++) {
    a2=XYZ(i,1) ; a3=a2;	  /* a2, a3 are min & max of XYZ(i) */
    for (iat=1;iat<=Num_atoms;iat++) {
	    aa=XYZ(i,iat) ; if (aa<a2) a2=aa ; if (aa>a3) a3=aa;
    } /* end for ( iat */
    abox0[i-1]=a2 ; abox1[i-1]=a3;
} /* end for  i */

for (i=0;i<3;i++) {
    abox0[i]=abox0[i] - aboxsz*0.5;
    Num_cubes[i]=(abox1[i]-abox0[i])/aboxsz;
    Num_cubes[i]++;
} /* end for ( i */
Number_of_planes=Num_cubes[0];
Number_of_lines=Num_cubes[0]*Num_cubes[1];

#define NBOX_Y Num_cubes[1]

xzero=abox0[0] ; yzero=abox0[1] ; zzero=abox0[2];

if (new) {
    fprintf(stderr,"HASH: lowest, highest coord, no of boxes\n");
    for (i=0;i<3;i++)
    fprintf(stderr,"HASH: %10.2f %10.2f %5d\n",
	abox0[i],abox1[i],Num_cubes[i]);
}

/*  create list of pointers and box indices */
/*  for each atom  Z_INDEX=z-box */
/*  NEXT=atom in same line with next higher (z) */
/*  for each box (line) FIRSTINLINE=index of atom with lowest (z) */
/*                      LOWZ_INLINE etc.=limits on z for a line */

/* fprintf(stderr,"HASH: num_lines_old,number_of_lines: %d %d\n",
    num_lines_old,number_of_lines); */

if (num_lines_old<Number_of_lines) {
    if (num_lines_old==0)
	AtomHash=(A_Hash *) malloc(Num_atoms*sizeof(A_Hash));
    else free(LineHash);

    LineHash=(L_Hash *) malloc(Number_of_lines*sizeof(L_Hash));

    num_lines_old=Number_of_lines;
}

/* initialize pointer arrays used for cubing */
for (ix=0;ix<Num_cubes[0];ix++) {
    for (iy=0;iy<Num_cubes[1];iy++) {
    FIRSTINLINE(ix,iy)= LOWZ_INLINE(ix,iy)= HIZ_INLINE(ix,iy)= -1;
}}
for (i=0;i<Num_atoms;i++) { Z_INDEX(i)= -1;  NEXT(i)= -1; }

x_pt=Xyz;
for (iat=0;iat<Num_atoms;iat++) {
    ix=(*(x_pt++)-xzero)/aboxsz;
    iy=(*(x_pt++)-yzero)/aboxsz;
    iz=(*(x_pt++)-zzero)/aboxsz;
    if (ix<0|ix>=Num_cubes[0] || iy<0|iy>=Num_cubes[1]) {
	fprintf(stderr,"HASH: ix or iy o.o.range= %4d %4d\n",ix,iy);
	exit(1);
    }
    /* Save box ind iz  "chain" */
    Z_INDEX(iat)=iz;
    where=FIRSTINLINE(ix,iy);
    lastwhere= -1;
    if (where== -1) {/* first box ever in this line */
        FIRSTINLINE(ix,iy)=iat;
        LOWZ_INLINE(ix,iy)=iz ; HIZ_INLINE(ix,iy)=iz ;
        continue;
    }
    while (TRUE) {/* find chains */
        if (iz<LOWZ_INLINE(ix,iy)) LOWZ_INLINE(ix,iy)=iz;
        if (iz>HIZ_INLINE(ix,iy)) HIZ_INLINE(ix,iy)=iz;
        bigger=FALSE ; equal=FALSE;
        if (iz>Z_INDEX(where)) bigger=TRUE;
        if ( bigger ) { lastwhere=where; /* keep looking */
            where=NEXT(where);
            if (where== -1) {
                NEXT(lastwhere)=iat ; break; } /* last in line */
        } /* end if bigger */
        else {/* insert before */
            NEXT(iat)=where;
            if (lastwhere== -1) FIRSTINLINE(ix,iy)=iat;
            else NEXT(lastwhere)=iat;
            break;
        } /* end insert before */
    } /* end repeat */
} /* end for ( iat */

if (new) {
    fprintf(stderr,"HASH: finished OK\n");
    /* This printout is too much when using LINES
    fprintf(stderr,"HASH: linefirst, hi & lowy, hi & loz=\n");
    for (ix=0;ix<nbox[0];ix++) {
	for (iy=0;iy<nbox[1];iy++) {
	    fprintf(stderr,"HASH: %3d %3d %5d%5d%5d\n",
	    ix,iy,FIRSTINLINE(ix,iy),LOWZ_INLINE(ix,iy),HIZ_INLINE(ix,iy));
    }}
    */
    new=FALSE;

#ifdef DEBUG
fprintf(stderr,"HASH: division according to lines:\n");
    for (i=0;i<Number_of_lines;i++) {
	fprintf(stderr,"line %d:",i);
	iat=FIRSTINLINE(i);
	while (iat>=0) {
	    fprintf(stderr," %d",iat); iat=NEXT(iat); }
	fprintf(stderr,"\n");
    }
/*    fprintf(stderr,"HASH: pointers by atom:\n");
*    for (i=0;i<nat;i++) { iat=i+1;
*	fprintf(stderr,"%3d %3d %3d %3d\n",iat,NEXT(i),Z_INDEX(i));
*    } */
#endif
}
} /* end HashXyz */

/* * * * * * * * * * * * * * * * * * * * * * * * * 
 *    Extracted from S*I*G*M*A - 1996 program
 * * * * * * * * * * * * * * * * * * * * * * * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * *
 *  ScanHash(xcop) 
 *    Use the hashing representation made by HashXyz to find 
 *    atoms within cutoff distance of xcop[3]
 * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "pair.H"

/* defined in header file "pair.H":
int BOXSIZE; set in input
int Num_box_cutoff,Num_cubes[],Number_of_planes; set by HashXyz()
*/

void ScanHash(REAL *xcop, int *pairs, int *npair, int Num_box_cutoff) 
{
int ix,iy,ix1,ix2,iy1,iy2,iz,iz1,iz2,iz3;
int where;
int notfirst;
int jat;
REAL dis2;

/* DETERMINE THE EXTENT OF HASH-SPACE WHERE NAYBORS MAY BE FOUND */
ix=((xcop[0]-xzero)/BOXSIZE); ix1=ix-Num_box_cutoff ; ix2=ix+Num_box_cutoff;
if (ix1>Num_cubes[0] || ix2<0) return;
iy=((xcop[1]-yzero)/BOXSIZE); iy1=iy-Num_box_cutoff ; iy2=iy+Num_box_cutoff;
if (iy1>Num_cubes[1] || iy2<0) return;
iz=((xcop[2]-zzero)/BOXSIZE); iz1=iz-Num_box_cutoff ; iz2=iz+Num_box_cutoff;
if (iz1>Num_cubes[2] || iz2<0) return;
if (ix1<0) ix1=0 ; if (ix2>=Number_of_planes) ix2=Number_of_planes-1;
if (iy1<0) iy1=0 ; if (iy2>=Num_cubes[1]) iy2=Num_cubes[1]-1;

/* ######### LOOP OVER LINE parallel to Z-axis, at ix,iy */
for (ix=ix1;ix<=ix2;ix++) {
for (iy=iy1;iy<=iy2;iy++) {
  where=FIRSTINLINE(ix,iy) ; if (where==-1) continue ; /* no points in line */
  if (iz1>HIZ_INLINE(ix,iy)) continue; if (iz2<LOWZ_INLINE(ix,iy)) continue;
  notfirst=FALSE;

/* CHAIN THRU THE ATOMS IN THIS LINE UNTIL ATOM-Z>IZ2 OR END OF LINE */
    while(TRUE) {
	if (notfirst) where=NEXT(where);
	else notfirst=TRUE;
	if (where==-1) break; /* done with line because no more points */
	iz3=Z_INDEX(where);
	/* determine if iz1 <=z<=iz2 */
	if (iz1>iz3) continue; if (iz2<iz3) break;

	/* ######## USE THE FOUND ATOM  ########## */
	jat=where+1;

	/* add to the pairlist */
	pairs[(*npair)++] = jat;

	/* For example: calculate the distance
	dis2= (xcop[0]-Xyz[3*jat])*(xcop[0]-Xyz[3*jat])+
	(xcop[1]-Xyz[3*jat+1])*(xcop[1]-Xyz[3*jat+1])+
	(xcop[2]-Xyz[3*jat+2])*(xcop[2]-Xyz[3*jat+2]);
	*/

    } /* end loop over atoms in line*/
}} /* end loop over line*/

} /* end ScanHash() */
