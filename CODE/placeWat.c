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
 *   placeWat program
 *
 *  Optimizes the water positions using the cavity coordinates as starting points.
 *
 *  Input:  
 *	argv[1] - protein structure
 *	argv[2] - starting oxygen placements
 *	argv[3] - minimization algorithm (optional - translate, rotate, or both - both is default)
 *  Output:  
 *	stdout  - minimized water placements
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

#define MIN_ITER_MAX 10       /* Max. # of iterations during BOTH min. algorithm */

extern void readPDB10(char *, int *, PDB **);
extern void initVDW (char *,int,PDB *);

extern REAL WaterEnergy   (REAL *);
extern REAL ElectEnergy   (REAL *);
extern REAL TransMin    (REAL *,REAL *);
extern REAL RotMin      (REAL *,REAL *);
extern REAL BothMin     (REAL *,REAL *);
extern REAL ScalarProduct3 ();
extern void CrossProduct(REAL v1[3], REAL v2[3], REAL result[3]);
extern int ConjGradMinimize(REAL *vars, int n_vars, REAL tol, int *iter, REAL *f_ret, REAL (*func)(REAL *, REAL *));
extern void OrientMolecule(REAL *angles);
extern void TranslateMolecule(REAL *xyz_new);
extern int PairList(REAL *ohh, int *NBPairs, int Num_atoms, PDB *pro);
extern void Align(int n, REAL *coords, int atom_idx, REAL *vec1, REAL *vec2);
extern void TurnMol(REAL angle, REAL *axis, REAL *coords, int start_idx, int num_atoms);

PDB *pro;                       /* protein coordinates */
int NBPairs[MAXPAIRS];          /* pair-list */
int NumberNBPairs;              /* size of pair-list */
REAL water[9];                  /* water coordinates */

void main(int argc, char *argv[])
{

char algorithm[10];             /* minimization algorithm (both,translate,or rotate) */
PDB *watO;                      /* input water coordinates */
REAL e_grad[9];                 /* electrostatic gradient for the water's oxygen */
REAL grad[3];                   /* potential energy gradient */
REAL MyGradient[6];             /* dummy gradient */
REAL RigidBodyVars[6];          /* variables used in minimization */
REAL fret;                      /* return value from minimizer */
REAL bisect[3];                 /* vector defining bisection of water */
REAL angle,angle_min;           /* angles used in coarse water placement */
REAL ener,ener_min;             /* potential energies */
REAL old_ener;                  /* energy used in iteration */
REAL d,d2;
FILE *fp;                       /* output file pointer */
int iter;                       /* number of conjugate gradient iterations */
int min_iter;                   /* number of minimization iterations */
int cgm_ok;			/* remember return code of minimizer */
int iwat;                       /* counter for input waters */
int numWat;                     /* number of waters */
int numPro;                     /* number of protein atoms */
int i,j;


if (argc < 3 || argc > 4)  {
    fprintf(stderr,"Correct usage: placeWat protein.pdb water.pdb min_algorithm (optional)\n\n");
    exit(1);
}

/* if min. algorithm not set (or incorrect) set default */
if (argc == 2)  {
    strcpy(algorithm,"both");
}
else  {
    strcpy(algorithm,argv[3]);
}
if ((EQU(algorithm,"both") != 1) && (EQU(algorithm,"trans") != 1) && (EQU(algorithm,"rot") != 1))  {
    fprintf(stderr,"Incorrect min. algorithm - use default (both)\n");
    strcpy(algorithm,"both");
}

/* read the two sets of coordinates */
readPDB10 (argv[1], &numPro, &pro);
readPDB10 (argv[2], &numWat, &watO);

/* read nonbonded parameters */
(void) initVDW (argv[1], numPro, pro);

/* clear the pairlist */
for(i=0;i<MAXPAIRS;i++) *(NBPairs+i)=0;  

for (iwat=0;iwat<numWat;iwat++)  {

/* Place oxygen on surface point */
    water[0] = watO[iwat].XX; water[1] = watO[iwat].YY; water[2] = watO[iwat].ZZ;
    NumberNBPairs = PairList(water,NBPairs,numPro,pro);
    
/* Initialize the coordinates */
    water[3] =  0.57735 + water[0]; 
    water[4] =  0.81650 + water[1]; 
    water[5] =            water[2]; 
    water[6] =  0.57735 + water[0]; 
    water[7] = -0.81650 + water[1]; 
    water[8] =            water[2]; 

/******************************************************************
 * Global placement of the water molecule:
 *   1. define electric field of oxygen
 *   2. place hydrogens in field
 *   3. rotate water around field vector by increment
 *   4. calculate electrostatic energy
 *   5. use placement at lowest electrostatic energy
 ******************************************************************/

/* ========================================================= */
/* Calculate electrostatic gradient for oxygen atom in water */
/* ========================================================= */
    ener_min = ElectEnergy(e_grad);
#ifdef DEBUG
    fprintf(stderr,"Initial Energy: %f\n",ener_min);
#endif

/* ================================================== */
/*          Align bisect and e_grad vectors           */
/* ================================================== */
    bisect[0] = -1; /* this sign change because e_grad has the wrong sign */
    bisect[1] =  0;
    bisect[2] =  0;
    Align(3,water,0,bisect,e_grad);

/* ==================================================
   Spin the water by 10 degree incrementss around the
   electrostatic gradient until elect. energy is min.
 * ================================================== */

    ener_min = WaterEnergy(e_grad);
#ifdef DEBUG
    fprintf(stderr,"Energy after field alignment: %f\n",ener_min);
#endif
    ener_min = 1.e10;
    angle_min = 0;
    angle = 10.*PI/180.;
    for (i=0;i<18;i++)  {
	ener = WaterEnergy(e_grad);
	if (ener < ener_min)  {
	    ener_min = ener;
	    angle_min = angle*i;
	}
	TurnMol(angle,bisect,water,0,3);
    }
    TurnMol(angle_min+PI,bisect,water,0,3);
#ifdef DEBUG
    ener_min = WaterEnergy(e_grad);
    fprintf(stderr,"Energy after spinning: %f\n",ener_min);
#endif

/******************************************************************
 * Refinement of water placement
 *   1. initialize coordinates to those from global placement
 *   2. initialize placement vector [1..6] to null vector
 *   3. calculate energies and gradients
 *   4. use minimizer to place new coordinates: rotation & translation
 *   5. iterate 3 & 4 until reach energy minimum
 ******************************************************************/

#ifndef FTOL
#define FTOL 1.0e-2
#endif

#ifdef DEBUG
	fprintf(stderr,"REMARK: input coordinates\n");
	fprintf(stderr,"ATOM%7i  OW  H2O%6i    %8.3f%8.3f%8.3f\n",3*iwat+1,iwat+1,water[0],water[1],water[2]);
	fprintf(stderr,"ATOM%7i  H1  H2O%6i    %8.3f%8.3f%8.3f\n",3*iwat+2,iwat+1,water[3],water[4],water[5]);
	fprintf(stderr,"ATOM%7i  H2  H2O%6i    %8.3f%8.3f%8.3f\n",3*iwat+3,iwat+1,water[6],water[7],water[8]);
#endif

/* ==================================================
   Minimize the energy using any of three options,
   translation alone, rotation alone or both,
   last one more than once.
 * ================================================== */

    for (j=0;j<3;j++) RigidBodyVars[j] = water[j];
    for (j=3;j<6;j++) RigidBodyVars[j] = 0.;
#ifdef DEBUG
    fprintf(stderr,"Energy before minimization = %f\n", BothMin(RigidBodyVars,MyGradient));
#endif

/***********************************************************************************************************
 *                              Minimization using translation                                             *
 ***********************************************************************************************************/
    if (EQU(algorithm,"trans"))  {
	cgm_ok = ConjGradMinimize(RigidBodyVars,3,FTOL,&iter,&fret,(*TransMin));
	TranslateMolecule(RigidBodyVars);
    }

/***********************************************************************************************************
 *                              Minimization using rotation                                                *
 ***********************************************************************************************************/
    else if (EQU(algorithm,"rot"))  {
	cgm_ok = ConjGradMinimize(RigidBodyVars+3,3,FTOL,&iter,&fret,(*RotMin));
	OrientMolecule(RigidBodyVars+3);
    }

/***********************************************************************************************************
 *                       Minimization using translation and rotation                                       *
 ***********************************************************************************************************/
    else if (EQU(algorithm,"both"))  {
	min_iter = 0;
	old_ener=1.e10;
	while (TRUE) {
	    for (j=0;j<3;j++) RigidBodyVars[j] = water[j];
	    for (j=3;j<6;j++) RigidBodyVars[j] = 0.;
	    NumberNBPairs = PairList(water,NBPairs,numPro,pro);
	    cgm_ok = ConjGradMinimize(RigidBodyVars,6,FTOL,&iter,&fret,(*BothMin));
	    TranslateMolecule(RigidBodyVars);
	    OrientMolecule(RigidBodyVars+3);
	    if (fabs(old_ener-fret) < 0.1) break;
	    if (min_iter > MIN_ITER_MAX)  {
		if (fret < old_ener) {
		    fprintf(stderr,
			"\n*** Energy did not converge in repeated minimization due to changes in the pairlist\n");
		    fprintf(stderr,
			"*** Test point #%d, after %d iterations energy difference = %f\n",
			iwat,MIN_ITER_MAX,fabs(old_ener-fret));
		    break;
		}
	    } 
	    old_ener=fret;
	    min_iter++;
	}
    }

    else {
	fprintf(stderr,"**** PlaceWat Error in selection of minimization algorithm\n");
	exit(1);
    }
#ifdef DEBUG
    fprintf(stderr,"Energy after minimization = %f\n",fret); 
    PrintMol(stderr);
#endif

    fprintf (stderr,".");

    /* test if the mol has moved "too far" */
    d = water[0] - watO[iwat].XX; d2  = d*d;
    d = water[1] - watO[iwat].YY; d2 += d*d;
    d = water[2] - watO[iwat].ZZ; d2 += d*d;
    if (d2 > MAXSHIFTSQ) continue;

    if (!cgm_ok) {
	fprintf (stderr,"\n*** Too many iterations in conjugate gradient minimization. ");
	fprintf (stderr,"Energy = %5.2f; test point #%d\n",fret,iwat);
    }
    fprintf(stdout,"HETATM%5d  OW  HOH%6d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
    watO[iwat].serial, watO[iwat].resSeq, water[0],water[1],water[2],
    WaterEnergy(MyGradient));
    fprintf(stdout,"HETATM%5d  H1  HOH%6d    %8.3f%8.3f%8.3f\n",
    watO[iwat].serial, watO[iwat].resSeq, water[3],water[4],water[5]);
    fprintf(stdout,"HETATM%5d  H2  HOH%6d    %8.3f%8.3f%8.3f\n",
    watO[iwat].serial, watO[iwat].resSeq, water[6],water[7],water[8]);

/*
    fprintf(stdout,"HETATM%5d  OW  HOH%6d    %8.3f%8.3f%8.3f\n",
    watO[iwat].serial, watO[iwat].resSeq, water[0],water[1],water[2]);
    fprintf(stdout,"HETATM%5d  H1  HOH%6d    %8.3f%8.3f%8.3f\n",
    watO[iwat].serial, watO[iwat].resSeq, water[3],water[4],water[5]);
    fprintf(stdout,"HETATM%5d  H2  HOH%6d    %8.3f%8.3f%8.3f\n",
    watO[iwat].serial, watO[iwat].resSeq, water[6],water[7],water[8]);
*/

   } /* end of loop over the input surface positions */

   exit (0);
} /* end of main() */
/* ********************************************************************** */


/******************************************************************
 *                                                                *
 *                          SUBROUTINES                           *
 *                                                                *
 ******************************************************************/


/******************************************************************
 *
 * PairList() - Build list of protein atoms within 12 angstroms
 *                of the water's oxygen atom
 *
 ******************************************************************/

int PairList(REAL *ohh, int *NBPairs, int Num_atoms, PDB *pro)

{
int numPairs,i,jat;
REAL aw[3],dist;

/* make a list of protein atoms near the water */
    numPairs=0-1; 
    for(jat=0;jat<Num_atoms;jat++) {
	for(i=0;i<3;i++) {
	    aw[0] = ohh[0] - pro[jat].XX;
	    aw[1] = ohh[1] - pro[jat].YY;
	    aw[2] = ohh[2] - pro[jat].ZZ;
	}
	dist = aw[0]*aw[0]+aw[1]*aw[1]+aw[2]*aw[2];
	if (dist > CUT_OFF*CUT_OFF) continue;

	/* meaning: atoms within 12A, but not itself */
	if (numPairs>=MAXPAIRS) {
	    printf ("Error: Num of NBPair atoms greater than MAXPAIRS parameter (%d)\n",MAXPAIRS);
	    exit(1);
	}
    
	numPairs++;
	NBPairs[numPairs] = jat;
    }

    numPairs++;
    return numPairs;

}  /* end PairList()   */


/******************************************************************
 *
 * PrintMol - Print water coordinates in PDB format
 *            (used in debugging)
 *
 ******************************************************************/

void PrintMol(FILE *file)
{

if (file)  { 
    fprintf(stderr,"No FILE* specified in PrintMol()\n");
    return;
} 
    fprintf(file,"ATOM%7i  OW  H2O%6i    %8.3f%8.3f%8.3f\n",1,1,water[0],water[1],water[2]);
    fprintf(file,"ATOM%7i  H1  H2O%6i    %8.3f%8.3f%8.3f\n",2,1,water[3],water[4],water[5]);
    fprintf(file,"ATOM%7i  H2  H2O%6i    %8.3f%8.3f%8.3f\n",3,1,water[6],water[7],water[8]);
}


/******************************************************************
 *
 * WaterEnergy() - Calculate total potential energy and its
 *                 gradient
 *
 *  Electrostatic:
 *      V = qi*qj/r
 *      dV/dx = -qi*qj*(xi-xj)/r3
 *
 *  VDW Energy:
 *      V = (1623*b)/r12 - (51*a)/r6 + qi*qj/r 
 *      dV/dx = [6*(a*51)/r6 - 12*(b*1623)/r12 - qi*qj/r] (xi-xj)/r2 
 *
 *  returns total energy in kcal/mol
 ******************************************************************/
   
REAL WaterEnergy(REAL grad[9])           /* total energy gradient */
{

int i,j,k,l;
REAL pcharge;           /* charge on the protein atom */
REAL dist[3];           /* distance (x,y,z) between protein and water atoms */
REAL wcharge[3];        /* charges on water atoms */
REAL r,r2;              /* distance between protein and water atoms */
REAL ener;              /* electrostatic energy */
REAL e_el;
REAL a,b,r6,r12,dudr2,e_lj6,e_lj12;
REAL sum_e_el,sum_e_lj;

wcharge[0] = C_OW*331;                    /* oxygen charge in water */
wcharge[1] = wcharge[2] =  C_H*331;       /* hydrogen charge in water */

for (i=0;i<9;i++)  grad[i] = 0.; 
sum_e_el = sum_e_lj = 0.;

/* Sum the contributions of each pair-list atom to the water's potential energy */
for (i=0;i<NumberNBPairs;i++) {
    j = NBPairs[i];
    pcharge = pro[j].LJ_c;
    a = pro[j].LJ_a; b = pro[j].LJ_b; 

    for (k=0;k<3;k++) {
        dist[0] = water[3*k]   - pro[j].XX;
        dist[1] = water[3*k+1] - pro[j].YY;
        dist[2] = water[3*k+2] - pro[j].ZZ;

        r2 = 1/(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2]);
	if (pcharge == 0.)  e_el = 0.;
	else  {
	    r = sqrt(r2);
	    e_el = pcharge*wcharge[k]*r;
	    sum_e_el += e_el;
	}	
	if (k == 0)  {
	    if ((a == 0)  && (b == 0))  {
		e_lj6 = e_lj12 = 0.;
	    }
	    else  {
		r6 = r2 * r2 * r2 ; r12 = r6 * r6;
		e_lj6 =   (A_OW * a * r6)/4.182;
		e_lj12 =  (B_OW * b * r12)/4.182;
		sum_e_lj += (e_lj12 - e_lj6);
	    }
	}
	else  e_lj6 = e_lj12 = 0.;

	dudr2 = (6. * e_lj6 - 12. * e_lj12 - e_el ) * r2;
	for (l=0;l<3;l++)  grad[l+3*k] += dudr2 * dist[l];
    } 
}

ener = sum_e_el + sum_e_lj;

#ifdef DEBUG
    fprintf(stderr,"ener: %f kcal/mol\n",ener);
    fprintf(stderr,"grad: %f %f %f %f %f %f %f %f %f\n",
	grad[0],grad[1],grad[2],grad[3],grad[4],grad[5],grad[6],grad[7],grad[8]);
#endif

return (ener);
}


/******************************************************************
 *
 * ElectEnergy() - Calculate electrostatic energy and its gradient
 *                 for the oxygen atom w/i water
 *                     V = qi*qj/r
 *                     dV/dx = -qi*qj*(xi-xj)/r3
 *
 *  -returns electrostatic energy in kcal/mol
 *  -Note: used only for global placement of the water 
 ******************************************************************/
  
REAL ElectEnergy(REAL e_grad[3])         /* elect. grad for oxygen atom */
{

int i,j,k,l;
REAL pcharge;           /* charge on the protein atom */
REAL dist[3];           /* distance (x,y,z) between protein and water atoms */
REAL wcharge[3];        /* charges on water atoms */
REAL r,r2;              /* distance between protein and water atoms */
REAL e_ener;            /* electrostatic energy */
REAL e_el;

wcharge[0] = C_OW*331;                    /* oxygen charge in water */
wcharge[1] = wcharge[2] =  C_H*331;       /* hydrogen charge in water */

for (i=0;i<9;i++)  e_grad[i] = 0.;
e_ener = 0.;

/* Sum the contributions of each NBPair atom to the water's elect. ener. */
for (i=0;i<NumberNBPairs;i++) {
    j = NBPairs[i];
    pcharge = pro[j].LJ_c;
    if (pcharge == 0)  continue; /* save time if charge == 0 */

    for (k=0;k<3;k++) {
        dist[0] = water[3*k]   - pro[j].XX;
        dist[1] = water[3*k+1] - pro[j].YY;
        dist[2] = water[3*k+2] - pro[j].ZZ;
        r2 = 1/(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2]);
        r = sqrt(r2);
        e_ener += (e_el = pcharge*wcharge[k]*r);
        for (l=0;l<3;l++) e_grad[3*k+l] -= e_el * dist[l] * r2;
    }
}

return (e_ener);
}


/******************************************************************
 *
 * RotGrad() - Calculate rotational gradient 
 *
 ******************************************************************/

int RotGrad(REAL grad[3], REAL rot_grad[9])       /* rotational gradient */

{
int i;
REAL torque1[3],torque2[3];
REAL spindle[3];
REAL r_oh[3];

for (i=0;i<3;i++)  r_oh[i] = (water[i+3] - water[i]);
CrossProduct(r_oh,grad+3,torque1);

for (i=0;i<3;i++)  r_oh[i] = (water[i+6] - water[i]);
CrossProduct(r_oh,grad+6,torque2);
for (i=0;i<3;i++) rot_grad[i] = -(torque1[i] + torque2[i]);

return 1;
}


/******************************************************************
 *
 * TransMin() - Minimize the water's potential energy using 
 *              translation 
 *
 ******************************************************************/

REAL TransMin(REAL *coor, REAL *RigidGrad)         /* energy grad. for min. algorithm */
{

REAL ener;               /* total potential energy for water */
REAL grad[9];       /* Cartesian grads for vdw and elect. energy */
REAL displ[3];
int i,iat;

REAL coor_save[9];
for (i=0;i<9;i++) coor_save[i] = water[i];

TranslateMolecule(coor);

ener = WaterEnergy(grad);

#ifdef DEBUG
fprintf(stderr,"TransMin Energy: %f\n",ener);
#endif

if (RigidGrad)  {   /* need grad so passed pointer otherwise grad==NULL */
    for (i=0;i<3;i++)  RigidGrad[i] = grad[i];
    for (i=0;i<3;i++) {
	    for (iat=1;iat<3;iat++) RigidGrad[i] += grad[3*iat+i];
	    RigidGrad[i] *= -1.;  /* ?????????????????????  */
    }
}

for (i=0;i<9;i++) water[i]=coor_save[i];

return (ener);
}


/******************************************************************
 *
 * TranslateMolecule - Move water molecule
 *
 ******************************************************************/

void TranslateMolecule(REAL *xyz_new)
{
REAL displ[3];
int i;

/* Minimizer uses oxygen coordinates for translational relaxation.
   Therefore, construct displacement vector for the entire water
   molecule from the old and new coordinates  */

for (i=0;i<3;i++)  {
    displ[i] = xyz_new[i] - water[i];
    water[i]   += displ[i];
    water[i+3] += displ[i];
    water[i+6] += displ[i];
}
} /* end TranslateMolecule() */


/******************************************************************
 *
 * RotMin() - Minimize the water's potential energy using
 *            rotation
 *
 ******************************************************************/

REAL RotMin(REAL *variables, REAL *rot_grad)         /* energy grad. for min. algorithm */
{

REAL ener;               /* total potential energy for water */
REAL vgrad[3],egrad[9];  /* grads for vdw and elect. energy */
REAL grad[9]; 
int i;

REAL coor_save[9];
for (i=0;i<9;i++) coor_save[i] = water[i];

/* rotate water around oxygen position */
OrientMolecule(variables);

ener = WaterEnergy(egrad);

#ifdef DEBUG
fprintf(stderr,"RotMin Energy: %f\n",ener);
#endif

if (rot_grad)  /* need grad so passed pointer otherwise grad==NULL */
    RotGrad(egrad,rot_grad);

for (i=0;i<9;i++) water[i] = coor_save[i];

return (ener);
}


/******************************************************************
 *
 * OrientMolecule - Rotate molecule around each coordinate axis
 *
 ******************************************************************/
void OrientMolecule(REAL *angles)
{
REAL spindle[3];

/* Rotate around x */
    spindle[0] = 1; spindle[1] = 0; spindle[2] = 0;
    TurnMol (angles[0], spindle, water, 0, 3);

/* Rotate around y */
    spindle[0] = 0; spindle[1] = 1; spindle[2] = 0;
    TurnMol (angles[1], spindle, water, 0, 3);

/* Rotate around z */
    spindle[0] = 0; spindle[1] = 0; spindle[2] = 1;
    TurnMol (angles[2], spindle, water, 0, 3);

} /* end OrientMolecule() */


/******************************************************************
 *
 * BothMin() - Minimize the water's potential energy using
 *             rotation and translation
 *
 ******************************************************************/

REAL BothMin(REAL *variables, REAL *RigidGrad)      /* gradient for rotation and translation */
{

REAL grad[9];        /* cartesian energy gradient */
REAL coor_save[9];   /* coordinates before rot. & trans. */
REAL ener;           /* total energy */
int i,iat;

for (i=0;i<9;i++) coor_save[i] = water[i];

TranslateMolecule(variables);
OrientMolecule(variables+3);

ener = WaterEnergy(grad);

#ifdef DEBUG
ifprintf(stderr,"BothMin Energy: %f\n",ener);
#endif


if (RigidGrad)  {   /* need grad so passed pointer otherwise grad==NULL */
    for (i=0;i<3;i++)  RigidGrad[i] = grad[i];
    for (i=0;i<3;i++) {
	for (iat=1;iat<3;iat++) RigidGrad[i] += grad[3*iat+i];
	RigidGrad[i] *= -1.;  /* ?????????????????????  */
    }
    
    RotGrad(grad,RigidGrad+3); /* the rotations are variables 3,4,5 */
}

/* Replace old coordinates so that water does not change yet */
for (i=0;i<9;i++) water[i] = coor_save[i];

return (ener);

} /* end BothMin() */
/* =================================================================== */
