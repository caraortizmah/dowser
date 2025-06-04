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
 *   refineWat program
 *
 *  Refines the treasure water placements by minimizing the energy of the waters.
 *  Note, the waters now "see" each other.
 *
 * Refinement of water placements:
 *   1. minimize energy of all water simulataneously
 *   2. iterate 1 until energy converges or reach max. allowed iterations 
 *
 *  Input:
 *      argv[1] - protein structure
 *      argv[2] - minimized water placements
 *      argv[3] - minimization algorithm (optional - translate, rotate, or both - both is default)
 *  Output:
 *      stdout - refined water placements
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

#define MAXPAIRS 1000                    /* Max. # of water-protein pairs */
#define MIN_ENERGY_DIFF 0.005            /* cutoff (fraction of change over energy) to stop energy minimization */
#define MAX_WATER_ITER 20                /* Max. # of minimization cycles */
#define WATER(i,j,k) water[9*i+3*j+k]

extern void readPDB10(char *, int *, PDB **);
extern void initVDW (char *,int,PDB *);

void appendWat(PDB *wat);
void WaterPairList(int Num_atoms);
REAL WaterEnergy(REAL *grad);
void RotGrad(REAL *grad, REAL *rot_grad, int water_num);
REAL TransMin(REAL *coor, REAL *RigidGrad);
void TranslateMolecule(REAL *xyz_new);
REAL RotMin(REAL *variables, REAL *rot_grad);
void OrientMolecule(REAL *angles);
REAL BothMin(REAL *variables, REAL *RigidGrad);

extern int ConjGradMinimize(REAL *p, int n, REAL ftol, int *iter, REAL *fret, REAL (*function_ptr)(REAL *, REAL *));
extern void CrossProduct(REAL *a, REAL *b, REAL *v);
extern void TurnMol(REAL angle, REAL *spindle, REAL *Xyz, int pivot, int Num_atoms);
REAL ScalarProduct3 ();
void PrintMol();

PDB *pro;                       /* protein coordinates */
int **NBPairs;                  /* pair-list */
int *NumberNBPairs;             /* size of pair-list */
REAL *water;                    /* water coordinates */
int numWat;                     /* number of waters to be refined */
int numProAtoms;                /* number of protein atoms */

int main(int argc, char *argv[])
{

char algorithm[10];             /* min. algorithm */
int i,j,count;                  /* counters */
int iwat;                       /* counter for input waters */
int iter;                       /* number of conjugate gradient iterations */
int min_iter;                   /* number of minimization iterations */
int cgm_ok;			/* remember return code of minimizer */
int numWatAtoms;                /* number of water atoms (should be 3*number of waters) */
int wat_iter;                   /* minimizations rounds on the water refinement */
int done;                       /* stop minimization? */
int allocation;                 /* memory allocation */
PDB *wat0;                      /* input water coordinates */
PDB *temp;                      /* temp. holder for PDB structure */
PDB *save_water;                /* temp. structure for storing old water info. */
REAL fret;                      /* return value from minimizer */
REAL bisect[3];                 /* vector defining bisection of water */
REAL angle,angle_min;           /* angles used in coarse water placement */
REAL ener,ener_min;             /* potential energies */
REAL old_ener;                  /* energy used in iteration */
REAL d,d2;                      /* distance the water moved during minimization */
REAL energy[MAX_WATER_ITER];    /* array storing the energy values for all of the waters */
REAL energy_diff;               /* array storing the energy change during minimization */
REAL *RigidBodyVars;            /* variables used in minimization */


#ifndef FTOL
    #define FTOL 1.0e-2
#endif

if (argc < 3 || argc > 4)  {
    fprintf(stderr,"USAGE: refineWat protein.pdb water.pdb min_alogrithm (optional)\n\n");
    exit(1);
}
/* if min. algorithm not set (or incorrect) set default */
if (argc == 3)  {
    strcpy(algorithm,"both");
}
else  {
    strcpy(algorithm,argv[3]);
}
if ((EQU(algorithm,"both") != 1) && (EQU(algorithm,"trans") != 1) && (EQU(algorithm,"rot") != 1))  {
    fprintf(stderr,"Incorrect min. algorithm - use default (both)\n");
    strcpy(algorithm,"both");
}

/* read the protein structure and the dowser waters */
readPDB10 (argv[1], &numProAtoms, &temp);
readPDB10 (argv[2], &numWatAtoms, &wat0);

numWat = (int) numWatAtoms/3;

/* temp is allocated for just numProAtoms.  However, we want more
   places to copy the waters into.  */
pro = (PDB *) malloc ((numProAtoms + numWatAtoms) * sizeof(PDB));
for (i=0;i<numProAtoms;i++)  {
    pro[i] = temp[i];
}
free (temp);

allocation = numWat*sizeof(REAL);
save_water    = (PDB *)  malloc (3*numWat*sizeof(PDB));
water         = (REAL *) malloc (9*allocation);
RigidBodyVars = (REAL *) malloc (6*allocation);
NBPairs       = (int **) malloc (numWat*sizeof(int *));
NumberNBPairs = (int *)  malloc (numWat*sizeof(int));

for (i=0;i<numWat;i++)  {
    NBPairs[i] = (int *)  malloc ((3*numWat+numProAtoms)*sizeof(int));
    /*
    NumberNBPairs[i] = 1;
    */
}

/* read nonbonded parameters for the protein structure */
(void) initVDW (argv[1], numProAtoms, pro);

/* Copy water into protein for pair-list */
appendWat(wat0);

/* Read initial coordinates into save_water */
for (i=0;i<3*numWat;i++)  {
    save_water[i] = wat0[i];
    water[3*i  ] = save_water[i].XX;
    water[3*i+1] = save_water[i].YY;
    water[3*i+2] = save_water[i].ZZ;
}

for (wat_iter=0;wat_iter<MAX_WATER_ITER;wat_iter++)  {

/***********************************************************************************************************
 *                              Initialize minimizer variables.                                            *
 *********************************************************************************************************** 
 *   Set first 3*numWat variables equal to the oxygen coordinates (translational relaxation).  Then,       * 
 *   initialize 3*numWat to 6*numWat equal to zero (rotational relaxation).                                *
 ***********************************************************************************************************/

    count = 0;
    for (j=0;j<numWat;j++)  {
	if (!strncmp(save_water[3*j].atomType,"O",1))  {
	    RigidBodyVars[count  ] = save_water[3*j].XX;
	    RigidBodyVars[count+1] = save_water[3*j].YY;
	    RigidBodyVars[count+2] = save_water[3*j].ZZ;
	    count += 3;
	}
    }
    for (j=3*numWat;j<6*numWat;j++)  RigidBodyVars[j] = 0.0;
    CALL WaterPairList(numProAtoms+3*numWat);

/***********************************************************************************************************
 *                              Minimization using translation                                             *
 ***********************************************************************************************************/
    if (EQU(algorithm,"trans"))  {
	cgm_ok = ConjGradMinimize(RigidBodyVars,3*numWat,FTOL,&iter,&fret,(*TransMin));
	TranslateMolecule(RigidBodyVars);
    }

/***********************************************************************************************************
 *                              Minimization using rotation                                                *
 ***********************************************************************************************************/
    else if (EQU(algorithm,"rot"))  {
	cgm_ok = ConjGradMinimize(RigidBodyVars+3*numWat,3*numWat,FTOL,&iter,&fret,(*RotMin));
	OrientMolecule(RigidBodyVars+3*numWat);
    }

/***********************************************************************************************************
 *                       Minimization using translation and rotation                                       *
 ***********************************************************************************************************/
    else if (EQU(algorithm,"both"))  {
	cgm_ok = ConjGradMinimize(RigidBodyVars,6*numWat,FTOL,&iter,&fret,(*BothMin));
	TranslateMolecule(RigidBodyVars);
	OrientMolecule(RigidBodyVars+3*numWat);
    }
    else {
	fprintf(stderr,"*** PlaceWat Error in selection of minimization algorithm\n");
	exit(1);
    }

#ifdef OLD
    if (!cgm_ok) {
	fprintf (stderr,"\n*** Too many iterations in conjugate gradient minimization. ");
	fprintf (stderr,"Energy = %5.2f; test point #%d\n",fret,iwat);
    }
#endif
    /* copy new water coordinates into set of current waters */
    for (i=0;i<3*numWat;i++)  {
	save_water[i].XX = water[3*i+0];
	save_water[i].YY = water[3*i+1];
	save_water[i].ZZ = water[3*i+2];
	save_water[i].tempFactor = pro[numProAtoms+i].tempFactor;
    }
    /* copy waters back into protein for new pair list calculation */
    appendWat(save_water);

/* stop minimizing if the difference in energy for all waters is less than MIN_ENERGY_DIFF */
    energy[wat_iter] = fret;
#ifdef DEBUG
    fprintf(stderr,"ITER: %d  ENERGY: %f\n",wat_iter,energy[wat_iter]);
#endif
    if (wat_iter > 0)  {
	energy_diff = fabs(energy[wat_iter] - energy[wat_iter-1])/fabs(energy[wat_iter]);
	if (energy_diff < MIN_ENERGY_DIFF)  break;
    }

}  /* end of minimization stage */

/* Ouput final water placements */
for (iwat=0;iwat<numWat;iwat++)  {
/*
    printf("HETATM%5d  OW  HOH W%4d    %8.3f%8.3f%8.3f%6.2f\n",
	wat0[3*iwat].serial, wat0[3*iwat].resSeq, save_water[3*iwat].XX,save_water[3*iwat].YY,
	save_water[3*iwat].ZZ, energy[iwat]);
*/
    printf("HETATM%5d  OW  HOH W%4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
	wat0[3*iwat].serial, wat0[3*iwat].resSeq, save_water[3*iwat].XX,save_water[3*iwat].YY,
	save_water[3*iwat].ZZ,save_water[3*iwat].tempFactor);
    printf("HETATM%5d  H1  HOH W%4d    %8.3f%8.3f%8.3f  1.00\n",
	wat0[3*iwat+1].serial, wat0[3*iwat+1].resSeq, save_water[3*iwat+1].XX,save_water[3*iwat+1].YY,
	save_water[3*iwat+1].ZZ); 
    printf("HETATM%5d  H2  HOH W%4d    %8.3f%8.3f%8.3f  1.00\n",
	wat0[3*iwat+2].serial, wat0[3*iwat+2].resSeq, save_water[3*iwat+2].XX,save_water[3*iwat+2].YY,
	save_water[3*iwat+2].ZZ);
}

return (wat_iter+1);
} /* end of main() */



/******************************************************************
 *                                                                *
 *                          SUBROUTINES                           *
 *                                                                *
 ******************************************************************/


/******************************************************************
 *
 * WaterPairList() - Build list of atoms within CUT_OFF angstroms from 
 *                   the oxygen atom of the water molecules.  Skip
 *                   the two hydrogens bonded to the oxygen
 *
 ******************************************************************/

void WaterPairList(int Num_atoms)        /* number of atoms in pro structure */
{

int numPairs,i,iat,jat;
int count;
REAL dist[3],dist2;

/* make a list of protein atoms near the waters' oxygen */
numPairs=0-1; 
for (iat=0;iat<numWat;iat++)  {
    count = 0;
    for(jat=0;jat<Num_atoms;jat++) {
	/* skip self */
	if (pro[jat].serial == pro[3*iat+numProAtoms].serial)  continue;

	dist[0] = pro[3*iat+numProAtoms].XX - pro[jat].XX;
	dist[1] = pro[3*iat+numProAtoms].YY - pro[jat].YY;
	dist[2] = pro[3*iat+numProAtoms].ZZ - pro[jat].ZZ;
	dist2 = dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2];
	if (dist2 > CUT_OFF*CUT_OFF) continue;

	/* If central atom is part of water (should be) then check whether the hydrogens for the
	   water are present.  If they are, exclude them from the pair-list.  */

	else if (!strcmp(pro[3*iat+numProAtoms].resName,"HOH"))  {
	    if ((pro[jat].serial == pro[3*iat+numProAtoms+1].serial) && !strcmp(pro[3*iat+numProAtoms+1].atomType,"H"))  {
		continue;
	    }
	    else if ((pro[jat].serial == pro[3*iat+numProAtoms+2].serial) && !strcmp(pro[3*iat+numProAtoms+2].atomType,"H"))  {
		continue;
	    }
	} 
	NBPairs[iat][count] = jat;
	count++;
    }   /* end of each atom loop */
    NumberNBPairs[iat] = count;
}       /* end of each water loop */

return;
}  /* end WaterPairList()   */


/******************************************************************
 *
 * appendWat() - Append dowser waters to the protein structure
 *
 ******************************************************************/

void appendWat(PDB *wat)
{

int i;

for (i=0;i<numWat;i++)  {
    wat[3*i  ].LJ_a = A_OW; wat[3*i  ].LJ_b = B_OW; wat[3*i  ].LJ_c = C_OW;
    wat[3*i+1].LJ_a = A_H;  wat[3*i+1].LJ_b = B_H;  wat[3*i+1].LJ_c = C_H;
    wat[3*i+2].LJ_a = A_H;  wat[3*i+2].LJ_b = B_H;  wat[3*i+2].LJ_c = C_H;
}
for (i=0;i<3*numWat;i++)  pro[i+numProAtoms] = wat[i];

return;
}  /* end of appendWat() */


/******************************************************************
 *
 * WaterEnergy() - Calculate total potential energy and its
 *                 gradient for all of the waters
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
   
REAL WaterEnergy(REAL *grad)           /* total energy gradient */
{

int i,iwat,j,k,l;
REAL dist[3];           /* distance (x,y,z) between protein and water atoms */
REAL wcharge[3];        /* charges on water atoms */
REAL r,r2;              /* distance between protein and water atoms */
REAL ener;              /* electrostatic energy */
REAL e_el;
REAL a,b,c,r6,r12,dudr2,e_lj6,e_lj12;
REAL sum_e_el,sum_e_lj;

wcharge[0] = C_OW*331;                    /* oxygen charge in water */
wcharge[1] = wcharge[2] =  C_H*331;       /* hydrogen charge in water */

for (i=0;i<3*numWat;i++)  grad[i] = 0.; 
ener = sum_e_el = sum_e_lj = 0.;

/* Sum the contributions of each pair-list atom to the waters' potential energy */
for (iwat=0;iwat<numWat;iwat++)  {
    for (k=0;k<3;k++) {
	for (i=0;i<NumberNBPairs[iwat];i++) {
	    j = NBPairs[iwat][i];
	    dist[0] = WATER(iwat,k,0) - pro[j].XX;
	    dist[1] = WATER(iwat,k,1) - pro[j].YY;
	    dist[2] = WATER(iwat,k,2) - pro[j].ZZ;
	    r2 = 1/(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2]);

	    a = pro[j].LJ_a; b = pro[j].LJ_b; c = pro[j].LJ_c; 
	    if (c == 0.)  e_el = 0.;
	    else  {
		r = sqrt(r2);
		e_el = c*wcharge[k]*r;
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
	    for (l=0;l<3;l++)  grad[3*iwat+l] += dudr2 * dist[l];
	}     /* end of pair-list loop */
#ifdef DEBUG
    if (iwat==2)  fprintf(stderr,"iwat: %d, k: %d, j: %d (%s), c: %f,  LJ: %f  EL: %f\n",iwat,k,j,pro[j].atomName,c,sum_e_lj,sum_e_el);
#endif
    }         /* end of ow,h1,h2 loop */
    pro[numProAtoms+3*iwat].tempFactor = sum_e_el + sum_e_lj;
    ener += sum_e_el + sum_e_lj;
    sum_e_el = sum_e_lj = 0.;

}             /* end of water loop */

#ifdef DEBUG
    fprintf(stderr,"LJ: %f  EL: %f\n",sum_e_lj,sum_e_el);
#endif
#undef DEBUG

return (ener);
}


/******************************************************************
 *
 * RotGrad() - Calculate rotational potential gradient.  Note, that
 *             the torque is defined as the negative rotational 
 *             gradient of the potential energy.  Calculate the
 *             torque using the bond from the oxygen to each hydrogen 
 *             as rotational arms.
 *
 ******************************************************************/

void RotGrad(REAL *grad, REAL *rot_grad, int water_num)        /* which water */
{
int i;
REAL torque1[3],torque2[3];
REAL spindle[3];
REAL r_oh[3];

for (i=0;i<3;i++)  r_oh[i] = (water[9*water_num+i+3] - water[9*water_num+i]);
CrossProduct(r_oh,grad,torque1);

for (i=0;i<3;i++)  r_oh[i] = (water[9*water_num+i+6] - water[9*water_num+i]);
CrossProduct(r_oh,grad,torque2);

for (i=0;i<3;i++) rot_grad[i] = -(torque1[i] + torque2[i]);

return;
}


/******************************************************************
 *
 * TransMin() - Calculate the energy and its gradient for 
 *              translational relaxation.  The water positions are 
 *              moved then the energy and gradient are calculated.
 *              Then the old coordinates are replaced.  
 *
 ******************************************************************/

REAL TransMin(REAL *coor, REAL *RigidGrad)    /* energy grad. for min. algorithm */
{

REAL ener;          /* total potential energy for water */
REAL *grad;         /* Cartesian grads for vdw and elect. energy */
REAL *coor_save;
int i,j,k,iat;

grad      = (REAL *) malloc (3*numWat*sizeof(REAL));
coor_save = (REAL *) malloc (9*numWat*sizeof(REAL));

for (i=0;i<9*numWat;i++)  coor_save[i] = water[i];

TranslateMolecule(coor);
ener = WaterEnergy(grad);

if (RigidGrad)  {   /* need grad so passed pointer otherwise grad==NULL */
    for (i=0;i<3*numWat;i++)  RigidGrad[i] = grad[i];
}

for (i=0;i<9*numWat;i++)  water[i] = coor_save[i];

free (grad);
free (coor_save);
return (ener);
}


/******************************************************************
 *
 * TranslateMolecule - Move all the water molecules based on the
 *                     new oxygen coordinates.  During 
 *                     translational relaxation, only the oxygen 
 *                     coordinates are used. 
 *
 ******************************************************************/

void TranslateMolecule(REAL *xyz_new)
{
REAL displ[3];
int i,j,k;

/* Construct displacement vector for the entire water
   molecule from the old and new oxygen coordinates  */

for (i=0;i<numWat;i++)  {
    for (j=0;j<3;j++)  {
	displ[j] = xyz_new[3*i+j] - WATER(i,0,j);
	for (k=0;k<3;k++)  {
	    WATER(i,k,j) += displ[j];
	}    /* end of each ow,h1,h2 */
    }        /* end of each x,y,z */
}            /* end of each water */

return;
} /* end TranslateMolecule() */


/******************************************************************
 *
 * RotMin() - Minimize the water's potential energy using
 *            rotation
 * TransMin() - Calculate the energy and its gradient for
 *              translational relaxation.  The water positions are
 *              moved then the energy and gradient are calculated.
 *              Then the old coordinates are replaced.
 *
 ******************************************************************/

REAL RotMin(REAL *variables, REAL *rot_grad)          /* energy grad. for min. algorithm */
{

REAL ener;               /* total potential energy for water */
REAL *grad; 
REAL *coor_save;
int i;

grad      = (REAL *) malloc (3*numWat*sizeof(REAL));
coor_save = (REAL *) malloc (9*numWat*sizeof(REAL));

for (i=0;i<9*numWat;i++) coor_save[i] = water[i];

/* rotate each water around oxygen position */
OrientMolecule(variables);

ener = WaterEnergy(grad);

if (rot_grad)  {   /* need grad so passed pointer otherwise grad==NULL */
    for (i=0;i<numWat;i++)  {
	RotGrad(grad+3*i,rot_grad+3*i,i);
    }
}
for (i=0;i<9*numWat;i++) water[i] = coor_save[i];

free (grad);
free (coor_save);
return (ener);
}


/******************************************************************
 *
 * OrientMolecule - Rotate each molecule around each coordinate axis
 *                  by a set of angles.
 *
 ******************************************************************/
void OrientMolecule(REAL *angles)
{
REAL spindle[3];
int iwat;

for (iwat=0;iwat<numWat;iwat++)  {
    /* Rotate around x */
	spindle[0] = 1; spindle[1] = 0; spindle[2] = 0;
	TurnMol (angles[3*iwat], spindle, water+9*iwat, 0, 3);

    /* Rotate around y */
	spindle[0] = 0; spindle[1] = 1; spindle[2] = 0;
	TurnMol (angles[3*iwat+1], spindle, water+9*iwat, 0, 3);
	/*TurnMol (angles[3*iwat], spindle, water+9*iwat, 0, 3);*/ 

    /* Rotate around z */
	spindle[0] = 0; spindle[1] = 0; spindle[2] = 1;
	TurnMol (angles[3*iwat+2], spindle, water+9*iwat, 0, 3);
	/*TurnMol (angles[3*iwat], spindle, water+9*iwat, 0, 3);*/
}

return;
} /* end OrientMolecule() */


/******************************************************************
 *
 * BothMin() - Minimize the water's potential energy using
 *             rotation and translation
 *
 ******************************************************************/

REAL BothMin(REAL *variables, REAL *RigidGrad)   /* gradient for rotation and translation */
{

REAL *grad;        /* cartesian energy gradient */
REAL *coor_save;   /* coordinates before rot. & trans. */
REAL ener;         /* total energy */
int i,iat;

grad      = (REAL *) malloc (6*numWat*sizeof(REAL));
coor_save = (REAL *) malloc (9*numWat*sizeof(REAL));

for (i=0;i<9*numWat;i++)  coor_save[i] = water[i];

TranslateMolecule(variables);
OrientMolecule(variables+3*numWat);

ener = WaterEnergy(grad);

if (RigidGrad)  {   /* need grad so passed pointer otherwise grad==NULL */
    for (i=0;i<3*numWat;i++)  RigidGrad[i] = grad[i];
    for (i=0;i<numWat;i++)  {
	RotGrad(grad+3*i,RigidGrad+3*i+3*numWat,i);
    }
}

/* Replace old coordinates so that water does not change yet */
for (i=0;i<9*numWat;i++)  water[i] = coor_save[i];

free (grad);
free (coor_save);
return (ener);
} /* end BothMin() */
/* ============================================================== */
