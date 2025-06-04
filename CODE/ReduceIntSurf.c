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
 *   reduceintsurf.c
 *      Reduces the number of internal surface points
 *
 *   Input:
 *      argv[1] - pdb-fomatted surface points (many)
 *   Output:
 *      argv[2] - ditto (few)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"
#define SEPARATION 1

extern void readPDB10(char *filename, int *num_atoms, PDB **atoms_ptr);
extern void WritePDB(FILE *fp, int nAtoms, PDB *atoms);

int main(int argc, char *argv[])
{
  int	numAtoms, numAtoms2=0;
  PDB	*pdbP, *pdb2;
  FILE *fp;
  int	i, j, k, found;
  REAL *xyz;
  REAL aaa, bbb;
  float separation;

  if (argc < 3 || argc > 4)  {
      fprintf(stderr,"USAGE: ReduceIntSurf reduces the number of points in an input file given a minimum separation\n");
      fprintf(stderr,"USAGE: ReduceIntSurf 'many_surf.pdb' 'few_surf.pdb' [SEPARATION]\n");
      exit(1);
  }
  if (argc == 4) sscanf (argv[3],"%f",&separation); 
  else separation = SEPARATION;
  fprintf (stdout, "* Reduce the number of internal surface points before running PlaceWat\n");
  fprintf (stdout, "* Separation between internal surface points will be %f\n", separation);
  separation *= separation;

  readPDB10(argv[1], &numAtoms, &pdbP);

  fprintf(stdout,"* Number of surface points in  input = %d\n", numAtoms);
  for (i=0; i<numAtoms; i++) {
      pdbP[i].xyz[0] = pdbP[i].XX;
      pdbP[i].xyz[1] = pdbP[i].YY;
      pdbP[i].xyz[2] = pdbP[i].ZZ;
  }

  pdb2 = (PDB *) malloc (numAtoms * sizeof (PDB));

  for (i=0; i<numAtoms; i++) {
      xyz =  pdbP[i].xyz;
      found = FALSE;
      for (j=0; j<numAtoms2; j++) {
	  aaa=0.;
	  for (k=0; k<3; k++) {
	      bbb=xyz[k] - pdb2[j].xyz[k];
	      aaa += bbb*bbb;
	  }
	  if (aaa<separation) { found=TRUE; break; }
      }

      if (!found) {
	 pdb2[numAtoms2].XX = pdb2[numAtoms2].xyz[0] = xyz[0];
	 pdb2[numAtoms2].YY = pdb2[numAtoms2].xyz[1] = xyz[1];
	 pdb2[numAtoms2].ZZ = pdb2[numAtoms2].xyz[2] = xyz[2];
	 numAtoms2++;
      }
  }
  fprintf(stdout,"* Number of surface points in output = %d\n", numAtoms2);


  fp = fopen (argv[2],"w");
  (void) WritePDB (fp, numAtoms2, pdb2);

  return 0;
}  /* end main */
