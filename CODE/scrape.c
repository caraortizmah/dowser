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
 *   Scrape.c
 *      Reject any point on surface A whose center is within a given distance of
 *      any single point on another surface, B.
 *
 *   Input:
 *      argv[1] - pdb filename containing Connolly surface with small probe
 *      argv[2] - pdb filename containing Connolly surface with large probe
 *   Output:
 *      stdout - remaining waters in PDB format
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

#define TRUE 1
#define FALSE 0

extern void readPDB10 (char *pdb_file, int *numAtoms, PDB **Atoms);
extern void WritePDB (FILE *fp, int nAtoms, PDB *atoms);

void main(int argc, char *argv[])
{
int	numWat, numBig, numSoak;
int	remain;
PDB	*wat, *big, *soakWat;
int	jat, iat, junk;
REAL	aa, bb, cc;
FILE	*fbigR;
REAL bigR, density;
char filename [256];

    strcpy(filename,getenv("DOWSER"));
    strcat(filename,"/DATA/ms_largeR.param");
    fprintf(stderr,"filename = %s\n", filename);
    fbigR = fopen(filename,"r");

    fscanf (fbigR, "%f%f%d%d", &density, &bigR, &junk, &junk);
    fclose (fbigR);
    fprintf (stderr, "bigR = %f\n", bigR);

    /* read surfaces A and B */
    readPDB10 (argv[1], &numWat, &wat);
    readPDB10 (argv[2], &numBig, &big);

    /* allocate space of the remaining sites */
    soakWat = (PDB *) malloc (numWat * sizeof(PDB));

    numSoak = 0;
    for (jat=1; jat<numWat; jat++) {
	remain = TRUE;
	for (iat=0; iat<numBig; iat++) {
	    aa = wat[jat].XX - big[iat].XX; aa *=aa; 
	    bb = wat[jat].YY - big[iat].YY; bb *=bb; 
	    cc = wat[jat].ZZ - big[iat].ZZ; cc *=cc; 
	    if ((aa+bb+cc) < bigR*bigR) {
		remain = FALSE;
		break;
	    }
	}
	if (remain) {
	    soakWat[numSoak] = wat[jat];
	    numSoak++; 
	}
    }

    fprintf (stderr, "numSoak = %d\n", numSoak );
    WritePDB (stdout, numSoak, soakWat);
}   /* end main() "scrape" */
