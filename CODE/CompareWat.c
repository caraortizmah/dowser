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
 *   CompareWat.c
 *      Compare water positions fiound by dowser and crystallographiuc water positions
 *
 *   Input:
 *      argv[1] - pdb file of dowser water molecules (tempFactor=energy)
 *      argv[2] - pdb file of xtal water molecules
 *      argv[3] - pdb file of internal xtal water molecules
 *   Output:
 *      report on standard output
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

extern void readPDB10 (char *, int *, PDB **);
REAL Distance2(PDB wat1, PDB wat2);

int main(int argc, char *argv[])

{
PDB *dowserwat_all, *dowserwat_int, *xtalwat_all, *xtalwat_int;
int nAtoms_dowser_all;
int nAtoms_dowser_int;
int nAtoms_xtal_all;
int nAtoms_xtal_int;

int i,j;
REAL dmin,dd;
int jmin;
PDB pdb1,pdb2;

    if (argc < 5 || argc > 6)  {
	fprintf(stderr,"USAGE: CompareWat dowser-wat_all.pdb dowser-wat_int.pdb xtal-wat_all.pdb xtal-wat_int.pdb\n");
	exit(1);
    }
    readPDB10 (argv[1], &nAtoms_dowser_all, &dowserwat_all);
    readPDB10 (argv[2], &nAtoms_dowser_int, &dowserwat_int);
    readPDB10 (argv[3], &nAtoms_xtal_all, &xtalwat_all);
    readPDB10 (argv[4], &nAtoms_xtal_int, &xtalwat_int);

    if ( nAtoms_dowser_int > 0 ) {
    fprintf (stdout,"\n* Find nearest xtal water for each dowser water\n");
    /* find nearest xtal water for each dowser water */
    fprintf (stdout,
    "     Dowser water         energy  distance       nearest xtal water\n");
    for (i=0;i<nAtoms_dowser_int;i+=3) {
	pdb1 = dowserwat_int[i];
	dmin=1.e4;
	for (j=0;j<nAtoms_xtal_all;j++) {
	    if ((dd = Distance2(pdb1,xtalwat_all[j])) < dmin) {
		dmin=dd; jmin=j;
	    }
	}
	pdb2=xtalwat_all[jmin];
	fprintf (stdout,
	    "#%2d %6.2f %6.2f %6.2f  %6.2f %6.2f A  #%3d %6.2f %6.2f %6.2f\n",
	    i/3+1,pdb1.XX,pdb1.YY,pdb1.ZZ,pdb1.tempFactor,
	    sqrt(dmin),pdb2.resSeq,pdb2.XX,pdb2.YY,pdb2.ZZ);
    }
    }

    if ( nAtoms_xtal_int > 0 ) {
    fprintf (stdout,"\n* Find nearest dowser water for each internal xtal water\n");
    /* find nearest dowser water for each xtal water */
    fprintf (stdout,
    "   Internal xtal water     energy distance       nearest Dowser water\n");
    for (i=0;i<nAtoms_xtal_int;i+=3) {
	pdb1 = xtalwat_int[i];
	dmin=1.e4;
	for (j=0;j<nAtoms_dowser_all;j+=3) {
	    if ((dd = Distance2(pdb1,dowserwat_all[j])) < dmin) {
		dmin=dd; jmin=j;
	    }
	}
	pdb2=dowserwat_all[jmin];
	fprintf (stdout,
	"#%3d %6.2f %6.2f %6.2f  %6.2f %6.2f A  #%2d %6.2f %6.2f %6.2f\n",
	pdb1.resSeq,pdb1.XX,pdb1.YY,pdb1.ZZ,pdb1.tempFactor,
	sqrt(dmin),jmin/3+1,pdb2.XX,pdb2.YY,pdb2.ZZ);
    }
    }

}

/* subroutine Distance2 */
REAL Distance2 (PDB wat1, PDB wat2)
{
    REAL dd,ddd;

    dd = wat1.XX - wat2.XX; ddd = dd*dd; 
    dd = wat1.YY - wat2.YY; ddd += dd*dd; 
    dd = wat1.ZZ - wat2.ZZ; ddd += dd*dd; 

    return ddd;
}
/* END */
