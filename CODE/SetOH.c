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
 *   setOH program
 *
 *  Optimize the hydrogen orientation in sidechain alcohols (i.e. in serine, threonine,
 *  and tyrosine).  The input pdb must be output from reformatPDB().
 *
 *  Input:
 *      argv[1] - reformatted protein structure (output of reformatPDB)
 *  Output:
 *      argv[2] - new protein structure
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

#define OXYGEN "OA"
#define SULFUR "S"
#define HYDROGEN "H"

extern void readPDB10(char *, int *, PDB **);
extern void initVDW (char *,int,PDB *);

int AlcoholPairList(int central_atom, int Num_atoms);
REAL MinimizeOhEnergy(REAL A[], REAL B[], int atom_num, int numPro);
REAL AlcoholEnergy(int atom_num, REAL grad[3]);
void UpdatePro(REAL *xyz, int atom_num);

extern void WriteXYZ(FILE *outfile, int num, PDB *atoms);
extern void TurnMol(REAL angle, REAL *spindle, REAL *Xyz, int pivot, int Num_atoms);

PDB *pro;                       /* protein coordinates */
int NBPairs[MAXPAIRS];          /* pair-list */
int NumberNBPairs;              /* size of pair-list */
REAL *water;

void main(int argc, char *argv[])
{

REAL energy_change = 0.;
REAL A[3],B[3];                 /* coordinates defining spindle */
int numPro;                     /* number of protein atoms */
int num_of_twist = 0;             /* number of serines in protein */
int i,j,ipro;
int back, backback;
FILE *outfile;                  /* FILE * to output pdb filename */

if (argc != 3)  {
    fprintf(stderr,"correct usage: setOH input.pdb output.pdb\n");
    exit(1);
}

outfile = fopen(argv[2],"w");
if (!outfile)  {
    fprintf(stderr,"*** Unable to open %s for output - use stdout\n",argv[2]);
    outfile = stdout;
}

/* read the protein coordinates */
readPDB10 (argv[1], &numPro, &pro);

/* read nonbonded parameters */
(void) initVDW (argv[1], numPro, pro);

/* clear the pairlist */
for(i=0;i<MAXPAIRS;i++) *(NBPairs+i)=0;  

/* loop over all the protein atoms */
for (ipro=0;ipro<numPro;ipro++)  {

/*******************************************************************
 * Find all hydrogens bonded to oxygen and sulfur
 *   via : atomtype == HYDROGEN, backchain is type OXYGEN or SULFUR
 *******************************************************************/
    if (!EQUAL(pro[ipro].type,HYDROGEN)) continue;
    if ((back=pro[ipro].back) == 0) continue;
    back = ipro - back;
    if ( EQUAL(pro[back].type,OXYGEN) ||
	 EQUAL(pro[back].type,SULFUR)) {
	num_of_twist++;
	for (i=0;i<3;i++) B[i] = pro[back].xyz[i];
	if ((backback=pro[back].back) == 0) continue;
	backback = back - backback;
	for (i=0;i<3;i++) A[i] = pro[backback].xyz[i];
        energy_change += MinimizeOhEnergy(A,B,ipro,numPro);
    }

}  /* end of protein atom loop */

WriteXYZ(outfile,numPro,pro);
fprintf(stdout,"* SetOH: Adjust %d O- and S- hydrogens on serine, threonine, tyrosine and cystine\n",
		 num_of_twist);
fprintf(stdout,"* SetOH: total energy change = %.2f\n", energy_change);

} /* end of SetOH() */


/******************************************************************
 *
 * REAL MinimizeOhEnergy() - Rotate alcohol hydrogen to find potential
 *                      energy minimium.
 *
 ******************************************************************/

REAL MinimizeOhEnergy(REAL A[], REAL B[], int atom_num, int numPro)              /* number of protein atoms */
{

REAL grad[3];            /* potential energy gradient */
REAL xyz[6];             /* coordinates of hydrogen and one back atom */
REAL spindle[3];         /* axis of rotation */
REAL ener,ener_min;      /* energy values */
REAL angle,angle_min;    /* angles */
int i;
REAL initial_energy;


/* define axis of rotation using the two back atoms from the hydrogen */
for (i=0;i<3;i++)  {
    spindle[i] = (B[i] - A[i]);
    xyz[i] = B[i];
    xyz[i+3] = pro[atom_num].xyz[i];  
}

NumberNBPairs = AlcoholPairList(atom_num,numPro);
initial_energy = ener_min = AlcoholEnergy(atom_num,grad);

#ifdef DEBUG
    fprintf(stderr,"Energy before spinning: %f\n",ener_min);
#endif

#define INCREMENTS 60   /* number of angles searched */
angle_min = 0;
angle = 2*PI/INCREMENTS;
for (i=0;i<INCREMENTS;i++)  {
    ener = AlcoholEnergy(atom_num,grad);
    if (ener < ener_min)  {
	ener_min = ener;
	angle_min = angle*i;
    }
    TurnMol(angle,spindle,xyz,0,2);
    UpdatePro(xyz+3,atom_num);
}

TurnMol(angle_min,spindle,xyz,0,2);
UpdatePro(xyz+3,atom_num);

    ener_min = AlcoholEnergy(atom_num,grad);
#ifdef DEBUG
    fprintf(stderr,"Energy after spinning: %f\n\n",ener_min);
#endif
    return (ener_min - initial_energy);

}  /* end of MinimizeOhEnergy() */


/******************************************************************
 *
 * UpdatePro() - Update the coordinates of a specified atom
 *
 ******************************************************************/

void UpdatePro(REAL *xyz, int atom_num)      /* atom number in pro structure */
{
int i;

for (i=0;i<3; i++) pro[atom_num].xyz[i] = xyz[i];

return;
} 


/********************************************************************
 *
 * AlcoholEnergy() - Calculate total potential energy and its
 *                   gradient (actually only the electrostatic b/c 
 *                   the vdw parameters for polar hydrogens are zero)
 *                   for alcohol sidechains.
 *
 *  Electrostatic:
 *      V = qi*qj/r
 *      dV/dx = -qi*qj*(xi-xj)/r3
 *
 *  returns total energy in kcal/mol
 ********************************************************************/
   
REAL AlcoholEnergy(int atom_num, REAL grad[3])           /* potential energy grad for hydrogen (not used) */
{

int i,j,k;
REAL c,c0;              /* charge on the protein atoms */
REAL dist[3];           /* distance (x,y,z) between protein and hydrogen */
REAL r,r2;              /* other distances */
REAL ener;              /* electrostatic energy */
REAL e_el;

for (i=0;i<3;i++)  grad[i] = 0.;
ener = 0.;

/* convert to kcal/mole */
c0 = pro[atom_num].LJ_c*331;

/* Sum the contribution of the pair-list atoms to the energy */
for (i=0;i<NumberNBPairs;i++) {
    j = NBPairs[i];
    c = pro[j].LJ_c;
    if (c == 0)  continue; /* save time if charge == 0 */

    for (k=0;k<3;k++) dist[k] = pro[atom_num].xyz[k] - pro[j].xyz[k];

    r2 = 1/(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2]);
    r = sqrt(r2);
    ener += (e_el = c*c0*r);
    for (k=0;k<3;k++) grad[k] -= e_el * dist[k] * r2;
}

return (ener);
}


/******************************************************************
 *
 * AlcoholPairList() - Build list of atoms within CUT_OFF angstroms
 *                  from a given protein atom
 *
 ******************************************************************/

int AlcoholPairList(int central_atom, int Num_atoms)        /* number of atoms in pro structure */
{

int numPairs,i,jat;
REAL dist[3],dist2;

/* make a list of protein atoms near the water */
numPairs=0-1;
for(jat=0;jat<Num_atoms;jat++) {

/* avoid self */
    if (jat == central_atom)  continue;
    /* Skip bonded atoms */
    else if (pro[jat].resSeq == pro[central_atom].resSeq)  {
	if      (!strcmp(pro[jat].resName,"SER"))  {
	    if (!strcmp(pro[jat].atomName,"OG"))  {
		continue;
	    }
	}
	else if (!strcmp(pro[jat].resName,"THR"))  {
	    if (!strcmp(pro[jat].atomName,"OG1"))  {
		continue;
	    }
	}
	else if (!strcmp(pro[jat].resName,"TYR"))  {
	    if (!strcmp(pro[jat].atomName,"OH"))  {
		continue;
	    }
	}
    }

    for (i=0;i<3;i++) dist[i] = pro[central_atom].xyz[i] - pro[jat].xyz[i];

    dist2 = dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2];
    if (dist2 > CUT_OFF*CUT_OFF) continue;

    numPairs++;

    if (numPairs>=MAXPAIRS) {
        printf ("Error: Num of NBPair atoms greater than MAXPAIRS parameter (%d)\n",MAXPAIRS);
        exit(1);
    }

    NBPairs[numPairs] = jat;
}

numPairs++;
return (numPairs);

}  /* end AlcoholPairList()   */
