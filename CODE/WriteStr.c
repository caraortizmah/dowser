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
 *
 *  WritePDB(): Write the protein str. in PDB format
 *  WriteXYZ(): Same as WritePDB but also outputs LJ parameters and charges
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

void WriteXYZ(FILE *outfile, int num, PDB *atoms);
void WritePDB(FILE *fp, int nAtoms, PDB *atoms);
void space(FILE *fp, int nSpace);

/***********************************************************************************
 *  WriteXYZ() : write the new pdb file with LJ parameters and charges
 ***********************************************************************************/

void WriteXYZ(FILE *outfile, int num, PDB *atoms)
{

int i;
PDB *anatom;
char aname[5],achar;

    for (i=0;i<num;i++) {
        anatom = atoms +i;
/*        if (!anatom->ftNote) continue;*/
        strcpy (aname,anatom->atomName);
        if (strlen(aname) == 4) {
            achar = aname[0];
            aname[0] = aname[1];
            aname[1] = aname[2];
            aname[2] = aname[3];
            aname[3] = NULLCHAR;
        }
        else achar = ' ';
        fprintf(outfile,
            "ATOM  %5d %c%-4s%-3s %5d   %8.3f%8.3f%8.3f  %4.2f%6.2f %7.3f %7.2f %7.1f\n",
            i+1, achar, aname, anatom->resName,
            anatom->resSeq, anatom->XX, anatom->YY, anatom->ZZ,
            anatom->occupancy, anatom->tempFactor,
            anatom->LJ_c, anatom->LJ_a, anatom->LJ_b);
    }
}


/***********************************************************************************
 *  WritePDB(): write the protein str. in PDB format w/occupanacy and beta values
 ***********************************************************************************/

void WritePDB (FILE *fp, int nAtoms, PDB *atoms)
{

  int i, nSpace;
  for (i=0; i<nAtoms; i++)  {
    nSpace = 6 - strlen (atoms[i].recdName);
    fprintf (fp, "%s", atoms[i].recdName);  space (fp,nSpace);

    fprintf (fp, "%5d", atoms[i].serial);

    space (fp,1);

    nSpace = 2 - strlen (atoms[i].atomType); space (fp,nSpace);
    fprintf (fp, "%s", atoms[i].atomType);

    nSpace = 2 - strlen (atoms[i].atomLoc);
    fprintf (fp, "%s", atoms[i].atomLoc);  space (fp,nSpace);

    space (fp,1);

    nSpace = 3 - strlen (atoms[i].resName);
    space (fp,nSpace);  fprintf (fp, "%s", atoms[i].resName);

    space (fp,1);

    nSpace = 1 - strlen (atoms[i].chainID); space (fp,nSpace);
    fprintf (fp, "%s", atoms[i].chainID);

    fprintf (fp, "%4d", atoms[i].resSeq);

    nSpace = 1 - strlen (atoms[i].iCode); space (fp,nSpace);
    fprintf (fp, "%s", atoms[i].iCode);

    space (fp,3);

    fprintf (fp, "%8.3f", atoms[i].XX);
    fprintf (fp, "%8.3f", atoms[i].YY);
    fprintf (fp, "%8.3f", atoms[i].ZZ);
    fprintf (fp, "%6.2f", atoms[i].occupancy);
    fprintf (fp, "%6.2f", atoms[i].tempFactor);

    space (fp,1);

    fprintf (fp, "%3d",atoms[i].ftNote);

    space (fp,2);

    nSpace = 4 - strlen (atoms[i].segID);
    fprintf (fp, "%s", atoms[i].segID);  space (fp,nSpace);

    nSpace = 2 - strlen (atoms[i].element);
    space (fp,nSpace);  fprintf (fp, "%s", atoms[i].element);

    nSpace = 2 - strlen (atoms[i].charge);
    fprintf (fp, "%s", atoms[i].charge);  space (fp,nSpace);

    fprintf (fp, "\n");
  }
  return;
}


/***********************************************************************************
 * space() - output a space to a FILE *
 ***********************************************************************************/

void space (FILE *fp, int nSpace)
{
  int j;
  for (j=0; j<nSpace; j++)
  fprintf (fp, " ");
}
