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
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *   chooser.c
 *      Choose a set of low energy non-overlaping O out of an overlaping set
 *   sparser.c
 *      Reduce the number of waters so that the oxygens of two waters are never closer 
 *      than 2.3 Angstroms.
 *
 *   Input:
 *      argv[1] - pdb file of water molecules (tempFactor=energy)
 *      argv[3] - selection criteria (distance, energy, or both (default))
 *   Output:
 *      argv[2] - pdb file containing remaining waters 
 *   Returns:
 *      number of waters remaining
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

#define	ALIVE	1
#define	DEAD	-1
#define MAX_ENERGY -10.0
#define	overlap	4.0
#define connect 2.3
#define contact connect*connect 

extern void readPDB10 (char *, int *, PDB **);
extern void WritePDB(FILE *fp, int nAtoms, PDB *atoms);
void chooser(int nWatAtom, PDB *wat, int *nAtoms_dense, PDB *dense);
void sparser(int nAtoms_dense, PDB *dense, int *nAtoms_sparse, PDB *sparse);
void SortByEnergy(int nWatAtom, PDB *water);

int main(int argc, char *argv[])

{
PDB *sparse, *wat, *dense;
FILE *fp;
char criteria[10];
int nWatAtom;
int nAtoms_sparse;
int nAtoms_dense;
int i,j;

    if (argc < 3 || argc > 4)  {
	fprintf(stderr,"USAGE: chooser input.pdb output.pdb selection_criteria (optional)\n");
	exit(0);
    }
    if (argc == 3)  {
       strcpy(criteria,"both");
    }
    else  {
       strcpy(criteria,argv[3]);
    }
    if (!strncmp(criteria,"dist",4) && !strncmp(criteria,"ener",4) && !strncmp(criteria,"both",4))  {
	fprintf(stderr,"Incorrect selection criteria - use default (both)\n");
	strcpy(criteria,"both");
    }

    readPDB10 (argv[1], &nWatAtom, &wat);

    sparse = (PDB *) malloc (nWatAtom * sizeof(PDB));

    /* Sort the waters by ascending energy */
    SortByEnergy(nWatAtom,wat);

    /* Use input criteria (or default) to "choose" waters */
    switch (criteria[0])  {

	case ('d'):     /* retain only non-overlapping set */
	    sparser(nWatAtom,wat,&nAtoms_sparse,sparse);
	    break;
	case ('e'):     /* retain only low energy set */
	    chooser(nWatAtom,wat,&nAtoms_sparse,sparse);
	    break;
	default:        /* retain both low energy and non-overlapping set */
	    dense = (PDB *) malloc (nWatAtom * sizeof(PDB));
	    chooser(nWatAtom,wat,&nAtoms_dense,dense);
	    sparser(nAtoms_dense,dense,&nAtoms_sparse,sparse);
	    break;
    }

    /* Sequentially order the Treasure waters */
    j = 0;
    for (i=0;i<nAtoms_sparse;i++)  {
	if (i%3 == 0) j++;
	sparse[i].serial = i+1;	          /* atom # */
	sparse[i].resSeq = j;             /* residue # */
	strcpy(sparse[i].chainID,"W");    /* chain ID */
    }
    fp = fopen(argv[2],"w");
    if (!fp)  {
	fprintf(stderr,"*** Unable to open file %s for writing - use stdout instead!\n",argv[2]);
	fp = stdout;
    }
    WritePDB (fp, nAtoms_sparse, sparse);
    return (nAtoms_sparse/3);

} /* end main() */


/*********************************************************************************
 * void SortByEnergy() - Sequentially sort all waters by ascending energy 
 *********************************************************************************/

void SortByEnergy(int nWatAtom, PDB *water)
{

PDB *temp;            /* temporary storage of water */
REAL low_ener;        /* lowest energy */
REAL next_ener;
int *index;           /* the ascending energy order */
int AllH;             /* logical for whether waters have hydrogens */
int nWat;             /* number of input waters */
int factor;
int i,j,k,count;


if ((water[0].resSeq==water[1].resSeq) && (water[0].resSeq==water[2].resSeq)) {
    AllH = TRUE;
    nWat = nWatAtom/3;
}
else {
    AllH = FALSE;
    nWat = nWatAtom;
}
temp  = (PDB *) malloc (nWatAtom*sizeof(PDB));
index = (int *) malloc (nWat*sizeof(PDB));

for (i=0; i<nWat; i++)  {
    index[i] = i;
    for (j=0; j<3; j++)  temp[3*i+j] = water[3*i+j];
}
if (AllH)  factor = 3;
else       factor = 1;

for (i=0; i<nWat; i++) {
    low_ener = temp[factor*index[i]].tempFactor;
    for (j=i+1;j<nWat;j++)  {
	next_ener = temp[factor*index[j]].tempFactor;
	if (next_ener < low_ener) {
	    low_ener = next_ener;
	    k = index[i]; index[i] = index[j]; index[j] = k;
	}
    }
}
for (i=0;i<nWat;i++)  {
    count = index[i];
    for (j=0;j<3;j++)  water[3*i+j] = temp[3*count+j];
}

free (temp);
return;
}


/*********************************************************************************
 * chooser() - retain waters with energies below -10 kcal/mol
 *********************************************************************************/

void chooser(int nWatAtom, PDB *wat, int *nAtoms_dense, PDB *dense)
{

int nWat;             /* number of input waters */
int AllH;             /* logical for whether waters have hydrogens */
int i,j,k;

if ((wat[0].resSeq==wat[1].resSeq) && (wat[0].resSeq==wat[2].resSeq)) {
    AllH = TRUE;
    nWat = nWatAtom/3;
}
else {
    AllH = FALSE;
    nWat = nWatAtom;
}

*(nAtoms_dense)=0;
for (i=0; i<nWat; i++) {
    if (AllH) {
	if (wat[3*i].tempFactor < MAX_ENERGY) {
	    dense[(*nAtoms_dense)++] = wat[3*i+0];
	    dense[(*nAtoms_dense)++] = wat[3*i+1];
	    dense[(*nAtoms_dense)++] = wat[3*i+2];
	}
    }
    else  {
	if (wat[i].tempFactor < MAX_ENERGY)  dense[(*nAtoms_dense)++] = wat[3*i];
    }
}

return;
} /* end chooser() */


/*********************************************************************************
 * sparser() - Retain a non-overlapping set of waters, preserve low energy waters
 *********************************************************************************/

void sparser(int nAtoms_dense, PDB *dense, int *nAtoms_sparse, PDB *sparse)
{

int jat,iat,i,j,k;
int near;              /* water is too close (or not) */
REAL aa, bb, cc;       /* distances */
PDB *anatom;           /* PDB pointer */

    if (nAtoms_dense == 0) {
	fprintf (stderr,"* No water molecules in the input of 'sparser'\n");
	exit (0);
    }

    /* The first water is always assumed to be in the non-conflicting set */
    for (i=0;i<3;i++) sparse[i] = dense[i]; 
    *nAtoms_sparse = 1; 
    for (jat=3; jat<nAtoms_dense; jat+=3) {
	near = FALSE;
	for (iat=0; iat<3* (*nAtoms_sparse); iat+=3) {
	    anatom = sparse + iat;
	    aa = dense[jat].XX - anatom->XX;
	    bb = dense[jat].YY - anatom->YY;
	    cc = dense[jat].ZZ - anatom->ZZ;
	    if ((aa*aa + bb*bb + cc*cc) < contact) { near = TRUE; break; }
	}
	if (near) continue;  /*then jat conflict with existing non-conflicting set*/
	/* Read the oxygen and two hydrogens into sparse */
	for (i=0;i<3;i++) { 
	    j=jat+i;
	    k=3* *nAtoms_sparse+i;
	    sparse[k] = dense[j]; 
	}
	(*nAtoms_sparse)++; 
    }
    (*nAtoms_sparse) *= 3;

return;
} /* end sparser */
