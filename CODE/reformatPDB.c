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
 *   reformatPDB.c
 *      Read in the original pdb file and a file with residue descriptions
 *      Eliminate nonpolar hydrogens, insert polar hydrogrens
 *      Choose proper chain termini (MUST DO)
 *      Detect residues that are not in the dictionary
 *      Read in additional data to help fix the problems (MUST DO)
 *
 *      Read in a file with geometry and connectivity, and compute missing xyz
 *
 *      Output a new pdb-formatted file
 *
 *   Input:
 *      -pdbin FILENAME - original pdb file 
 *      -atomtypes FILENAME or
 *          $DOWSER/DATA/atomdict.db - dictionary file with residue descriptions
 *      -atomparms FILENAME or
 *          $DOWSER/DATA/atomparms.db - dictionary file with LJ parameters 
 *   Output:
 *      -pdbout FILENAME - new PDB file
 *   Additional optional arguments
 *      -hetero = use ATOM and HETATM records in the pdb file
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

#define SSBONDSQ  5.0 /* if S-S distance^2 is below this, it's a bond */
#define NCBONDSQ  4.0 /* if N-C distance exceeds this, it's a chain break */
#define MAXTYPES 500
#define MAXMERGE 100

#define ERR_GENERAL 1
#define ERR_DOWSER_ENV_DICT 2
#define ERR_OPEN_ATOMDICT 3
#define ERR_DOWSER_ENV_PARAM 4
#define ERR_OPEN_ATOMPARMS 5
#define ERR_SCANDICTPDB 6
#define ERR_ATOM_COORDS 7
#define ERR_MAXTYPES_EXCEEDED 8
#define ERR_MAXMERGE_EXCEEDED 9

/* begin0,end0 = begin and end  in the set that remains when heteroatoms are removed 
   begin,end   = begin and end  in the set containing also the heteroatoms
   newbegin,newend = begin and end in ATOMTYPE */
typedef struct {
  char name[5],terminus[5];
  int begin0,end0,begin,end,newbegin,newend,type,numat,numextra;
  int chainend,index; /* 1 for N-term, 2 for C-term */
} RESIDUE;

extern void  readPDB10(char *pdbFileName, int *nAtomP, PDB **atomsP); /* Assuming from dowser.h, params named for clarity */
extern int   readDict (char *filename, int *numat, int *numrestyp, ATOMTYPE **atomtypes, RESTYPE **restypes);
extern int   BackOneBond(int from, PDB *atoms, ATOMTYPE *atomtypes);
extern char  *getenv(const char *name); /* Standard prototype */
extern int   ForwardOneBond(int from, PDB *atoms, ATOMTYPE *atomtypes, int nAtoms);
extern void  FindSSBonds (FILE *outFileP, PDB *atoms, RESIDUE *residues, int nRes, RESTYPE *restypes, int nRestype);
extern int   CysSG (PDB *atoms, int resnum, RESIDUE *residues, REAL *xs1);
extern int   Add1Atom(PDB *atoms, int backback_idx, int backatom_idx, int refatom_idx, int atom_to_place_idx, REAL bond, REAL angle, REAL dihedral); /* Updated to guessed ANSI C prototype */
extern void  FindChainBreaks(FILE *outFileP, PDB *atoms, RESIDUE *residues, int nRes, RESTYPE *restypes, ATOMTYPE *atomtypes);
extern int   ScanDictVsPDB(int nRes, RESIDUE *residues, int nRestype, RESTYPE *restypes);
extern int   MergeTerminus(FILE *outFileP, int countonly, PDB *atoms, RESIDUE *aresidue, int nRestype, RESTYPE *restypes, ATOMTYPE *atomtypes, ATOMTYPE *mergeatomtype, RESTYPE *mergerestype);
extern int   readParam (char *filename, int nAtoms, PDB *atoms);
extern void  WriteXYZ (FILE *outFileP, int num, PDB *atoms, RESIDUE *residues);
extern int   LocateInResidue(PDB *atoms, int resnum, RESIDUE *residues, char *name); /* Added extern */
extern int   AddAllAtoms (int nAtoms, PDB *atoms, ATOMTYPE *atomtypes, int final_pass); /* Added extern */

/* FILE *outfile; */ /* Removed global variable */

ATOMTYPE *atomtypes;

int main(int argc, char *argv[])
{
  FILE *outfile = NULL; /* outfile is now local to main */
  PDB	*newatoms, *atoms, *anatom;
  RESIDUE *residues, *aresidue;
  ATOMTYPE *anatomtype, /* *bnatomtype, */ *useatomtypes; /* bnatomtype removed, only used in #ifdef TESTM */
  RESTYPE *restypes, *arestype;
  ATOMTYPE mergeatomtype[MAXMERGE];
  RESTYPE mergerestype;
  /* int *type_key; */ /* Unused */
  int nType, nRestype;

  int i, j, k, ires, iat;
  int nAtom, nRes;
  int newnAtoms;
  int oldresSeq; char oldchainID; char oldinsert;
  /* char *oldtype; */ /* Unused */
  char *what; /* Used for reading atom names from anatomtype->name */
  int found; /* Used for searching atoms */
  char junk [10];
  int pdbin = 0, pdbout =0;
  int error=0;

  int noheteros = TRUE;
  char filetypes[256],fileparms[256];

  /* look for arguments that affect the way dowser is to be used */
  *filetypes = *fileparms = NULLCHAR;
  for (i=1; i<argc; i++) {
      if (EQUAL(argv[i],"-hetero")) noheteros = FALSE;
      if (EQUAL(argv[i],"-atomtypes")) {
          strncpy(filetypes,argv[i+1], sizeof(filetypes) - 1);
          filetypes[sizeof(filetypes)-1] = NULLCHAR;
      }
      if (EQUAL(argv[i],"-atomparms")) {
          strncpy(fileparms,argv[i+1], sizeof(fileparms) - 1);
          fileparms[sizeof(fileparms)-1] = NULLCHAR;
      }
      if (EQUAL(argv[i],"-pdbin")) pdbin=i+1;
      if (EQUAL(argv[i],"-pdbout")) pdbout=i+1;
  }

  if (pdbin*pdbout == 0) {
      fprintf(stderr,"ERROR: reformatPDB needs -pdbin FILENAME -pdbout FILENAME\n");
      return ERR_GENERAL;
  }

  if (noheteros) fprintf(stderr,"* DOWSER - default, no HETATM records will be used\n");
  else {
      fprintf(stderr,"* DOWSER - both ATOM and HETATM records will be used\n%s\n",
      		     "*          except atom 'O' in residue 'HOH'");
  }

  /* process the original PDB file */
  /* the readPDB routine allocates the needed storage for atoms[] */
  readPDB10(argv[pdbin], &nAtom, &atoms);

  if (!(outfile = fopen (argv[pdbout],"w"))) {
    fprintf (stderr,"REMARK ERROR: cannot open output pdb file %s\n",argv[pdbout]);
    return ERR_GENERAL;
  }


/*
ATOM     60  C   ASN    25      33.889  33.642  29.803  1.00 12.34           C  
*/
  /* count the number of residues in the ATOM records of the PB file, based on changes in
	 chain ID, residue sequence number, or insertion code. This is a first pass to determine nRes. */
  anatom = atoms;
  nRes = 0; oldresSeq=0; oldchainID='\0'; oldinsert=' ';
  for (iat=0; iat<nAtom; iat++) {

    if ( !EQUAL(anatom->recdName,"ATOM") &&
         !EQUAL(anatom->recdName,"HETATM")) { anatom++; continue; }
    if ( noheteros ) {
	if (EQUAL(anatom->recdName,"HETATM")) {
	    if (!EQUAL(anatom->resName,"HOH")) {
		fprintf(stderr,"* HETERO atom %s %s %s %d%s will be skipped\n",
		anatom->atomType, anatom->resName, anatom->chainID, anatom->resSeq, anatom->iCode);
	    }
	    anatom++; continue;
	}
    }
    else if (EQUAL(anatom->recdName,"HETATM") && EQUAL(anatom->resName,"HOH"))
	{ anatom++; continue; }

    if (oldresSeq != anatom->resSeq || oldchainID != anatom->chainID[0] || oldinsert != anatom->iCode[0])
	{ nRes++ ; oldresSeq = anatom->resSeq; oldchainID = anatom->chainID[0];
	  oldinsert = anatom->iCode[0]; }
    anatom++;
  }
  fprintf (outfile,"REMARK Number of residues = %d\n",nRes);

  /* allocate the structure describing the residues in the input PDB */
  residues = (RESIDUE *) malloc (nRes * sizeof(RESIDUE));

  /* fill the Residue array (residues[]) with information from the input PDB (atoms[]).
     This pass identifies residue boundaries and names. */
  anatom = atoms; aresidue = residues;
  nRes = 0; oldresSeq=0; oldchainID='\0';
  for (iat=0; iat<nAtom; iat++) {

    if ( !EQUAL(anatom->recdName,"ATOM") &&
         !EQUAL(anatom->recdName,"HETATM")) { anatom++; continue; }
    if ( noheteros ) {
	if (EQUAL(anatom->recdName,"HETATM")) { anatom++; continue; }}
    else if (EQUAL(anatom->recdName,"HETATM") && EQUAL(anatom->resName,"HOH"))
	{ anatom++; continue; }
	  

    /* new residue if residue number or chain-id changes */
    if (oldresSeq != anatom->resSeq || oldchainID != anatom->chainID[0] || oldinsert != anatom->iCode[0])
	{
	    nRes++ ; oldresSeq = anatom->resSeq; oldchainID = anatom->chainID[0];
	    oldinsert = anatom->iCode[0];
	    strncpy(aresidue->name, anatom->resName, sizeof(aresidue->name) - 1);
        aresidue->name[sizeof(aresidue->name)-1] = NULLCHAR;
	    aresidue->begin0 = aresidue->end0 = iat; /* init. end for 1-atom residues */
	    aresidue->chainend = aresidue->numextra = aresidue->numat = 0;
	    aresidue->index = nRes;
	    aresidue++;
	}
    else {
	(aresidue-1)->end0=iat;
    }
    anatom++;
  }
  fprintf(stdout,"* Input pdb file contains %d atoms\n",(int)(anatom-atoms));

#ifdef TEST
  fprintf (stdout,"Residues\n");
  for (i=0;i<nRes;i++) fprintf (stdout,"RES: %s %d %d %d %d\n",residues[i].name,
	  residues[i].begin, residues[i].end,
	  residues[i].begin0, residues[i].end0);
#endif

  /* read in the file describing the atoms with charges, connectivity, etc., by residue */
  if (readDict (filetypes, &nType, &nRestype, &atomtypes, &restypes) != 0) {
      /* Error message already printed by readDict */
      free(residues); /* Clean up allocated memory so far */
      /* atoms is allocated by readPDB10, needs to be handled if readPDB10 is also changed */
      if (outfile) fclose(outfile);
      return ERR_GENERAL; /* Or a more specific error from readDict if desired */
  }

  /* scan the residue types in the molecule vs. those in the dictionary */
  if (ScanDictVsPDB(nRes,residues,nRestype,restypes)) {
      fprintf(stderr, "ERROR: ScanDictVsPDB failed.\n");
      free(residues);
      free(atomtypes);
      free(restypes);
      if (outfile) fclose(outfile);
      return ERR_SCANDICTPDB;
  }

  /* Locate SS bridges and change the residue names to CSS */
  (void) FindSSBonds (outfile, atoms,residues,nRes,restypes,nRestype);

  /* Locate chain breaks and mark these in the residues[] by setting terminus */
  (void) FindChainBreaks(outfile, atoms,residues,nRes,restypes,atomtypes);

  /* Figure the extra atoms needed for the chain termini.
     This is done by calling MergeTerminus in 'countonly' mode (TRUE).
     It updates numextra in each RESIDUE struct. */
  for (ires=0;ires<nRes;ires++) {
      if (residues[ires].chainend && *(residues[ires].terminus)) { /* If it's a defined chain terminus */
	  if (MergeTerminus(outfile, TRUE,atoms,residues+ires,nRestype,restypes,atomtypes,
	      mergeatomtype,&mergerestype) != 0) {
              error = ERR_MAXMERGE_EXCEEDED;
              fprintf(stderr, "Error: MAXMERGE limit reached during atom counting for terminus.\n");
              /* Free memory before exiting */
              free(residues);
              free(atoms); /* Assuming atoms is always allocated by readPDB10 if program reaches here */
              if (atomtypes) free(atomtypes);
              if (restypes) free(restypes);
              if (outfile) fclose(outfile);
              return error;
          }
      }
  }

  /* Compute the total number of atoms in the new molecule (newnAtoms)
     by summing atoms from each residue type in the dictionary plus any extra atoms for termini.
     Also sets 'begin' and 'newend' indices for each residue in the context of the newatoms array. */
  newnAtoms = 0;
  for (ires=0; ires<nRes; ires++) {
	  aresidue = residues+ires;
	  aresidue->begin = newnAtoms; /* Start index of this residue in the newatoms array */
	  arestype = restypes + aresidue->type; /* Get the residue type from dictionary */
      /* Number of atoms for this residue = (atoms in dict definition) + (extra atoms for terminus) */
	  newnAtoms += arestype->end - arestype->begin + 1 + aresidue->numextra;
	  aresidue->newend = newnAtoms-1; /* End index of this residue in the newatoms array */
  }
  fprintf (outfile,"REMARK number of atoms in output molecule = %d\n",newnAtoms);
  newatoms = (PDB *) malloc (newnAtoms * sizeof (PDB)); /* Allocate memory for the new atom list */

  /* Fill the newatoms array.
     This loop iterates through each residue identified in the input PDB.
     For each residue, it iterates through the atoms defined in its (possibly modified for terminus) dictionary type.
     It tries to find a matching atom in the input PDB to copy coordinates and other info.
     If not found, it takes the atom name from the dictionary and marks it for coordinate calculation.
   */
  newnAtoms=0; /* Reset newnAtoms to use as an index for filling the newatoms array */
  for (ires=0; ires<nRes; ires++) {
      aresidue = residues + ires;
      arestype = restypes + aresidue->type; /* Default residue type */

      /* is there a terminus to add? If so, MergeTerminus is called in 'actual merge' mode (FALSE).
         This will populate mergeatomtype[] and mergerestype with the combined atom definitions.
         arestype and useatomtypes are then pointed to these merged definitions for this residue. */
      if (aresidue->chainend && *(aresidue->terminus)) {
	  if (MergeTerminus(outfile, FALSE,atoms,residues+ires,nRestype,restypes,atomtypes,
	      mergeatomtype,&mergerestype) != 0) {
              error = ERR_MAXMERGE_EXCEEDED;
              fprintf(stderr, "Error: MAXMERGE limit reached during atom merging for terminus.\n");
              /* Free memory before exiting - newatoms might not be fully populated or even allocated if this is the first error */
              free(residues);
              free(atoms);
              if (newatoms) free(newatoms); 
              if (atomtypes) free(atomtypes);
              if (restypes) free(restypes);
              if (outfile) fclose(outfile);
              return error;
          }
	  arestype=&mergerestype;
	  useatomtypes=mergeatomtype;

#ifdef TESTM
    anatomtype=mergeatomtype;
    for (i=0;i<mergerestype.numat; i++) {
	j=anatomtype->index;
	bnatomtype=atomtypes+j;
	fprintf(stdout,"%4d %-4s %-4s %-4s %-4s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
	  j,bnatomtype->name,bnatomtype->back,bnatomtype->forward,bnatomtype->type,
	  bnatomtype->charge,bnatomtype->bond,bnatomtype->angle,bnatomtype->dihedral);
	anatomtype ++;
    }
#endif
      }
      else useatomtypes=atomtypes;

      for (j=arestype->begin; j<=arestype->end; j++) {
	  anatomtype = useatomtypes + j; 
	  what = anatomtype->name; /* Name of the current atom from the dictionary definition */
	  found = 0; /* Flag to indicate if this dictionary atom is found in the input PDB residue */
          /* Search for the current dictionary atom ('what') within the original atoms of this residue */
	  for (k=aresidue->begin0; k<= aresidue->end0; k++) {
	      strncpy(junk,atoms[k].atomName, sizeof(junk) -1); 
          junk[sizeof(junk)-1] = NULLCHAR;
	      if (EQUAL(what,junk))  { 
		  found = 1; 
		  break;  /* Atom found in input PDB */
	      }
	  }
	  if (found) { /* Atom from dictionary was found in the input PDB */
	      newatoms[newnAtoms] = atoms[k]; /* Copy the PDB record (coords, B-factor, etc.) */
	      newatoms[newnAtoms].ftNote = TRUE; /* Mark that original coordinates are present */
	      /* fprintf(stdout,"FOUND %s in %d %s\n",what,i,aresidue->name); */
	  }
	  else { /* Atom from dictionary was NOT found in the input PDB (e.g., a hydrogen to be added) */
	      strncpy(newatoms[newnAtoms].atomName,what, sizeof(newatoms[newnAtoms].atomName) - 1);
          newatoms[newnAtoms].atomName[sizeof(newatoms[newnAtoms].atomName)-1] = NULLCHAR;
	      strncpy(newatoms[newnAtoms].resName,arestype->resname, sizeof(newatoms[newnAtoms].resName) - 1);
          newatoms[newnAtoms].resName[sizeof(newatoms[newnAtoms].resName)-1] = NULLCHAR;
	      newatoms[newnAtoms].ftNote = FALSE;
	      if (what[0]!='h' && what[0]!='H')
		  fprintf(outfile,"REMARK WARNING non-hydrogen NOT FOUND %s in %d %s\n",
		      what,ires+1,aresidue->name);
	  }
      /* Common properties for both found and new atoms, mostly from dictionary */
	  newatoms[newnAtoms].newResNo = ires; /* Set new residue number (0-indexed) */
	  newatoms[newnAtoms].LJ_c     = anatomtype->charge; /* Copy charge from dictionary */
	  newatoms[newnAtoms].key_dict = anatomtype->index;  /* Link to the master atomtype entry */
	  strncpy(newatoms[newnAtoms].atomType, anatomtype->type, sizeof(newatoms[newnAtoms].atomType) - 1);
      newatoms[newnAtoms].atomType[sizeof(newatoms[newnAtoms].atomType)-1] = NULLCHAR;
	  newnAtoms++;
      }

  } /* end fill new array of atoms */

#ifdef TESTX
  for (i=0;i<newnAtoms;i++) {
	  fprintf(stdout,"NEWATOM %-4s %-4s ",newatoms[i].resName,newatoms[i].atomName);
	  if (newatoms[i].ftNote) fprintf(stdout," OLD "); /* ftNote TRUE means coords from input PDB */
	  else fprintf(stdout," NEW "); /* ftNote FALSE means coords to be calculated */
	  j=newatoms[i].key_dict;
	  fprintf(stdout," %4d %-4s %-4s %-4s\n",
	      j,atomtypes[j].back,atomtypes[j].forward,atomtypes[j].type);
  }
#endif

    /* Iteratively add missing atom coordinates.
       AddAllAtoms is called multiple times because the calculation of one atom's coordinates
       might depend on another atom whose coordinates were calculated in a previous iteration.
       The loop continues until no new atoms can be placed or an error/stall is detected.
       The 'final_pass' argument (TRUE for the last call in an iteration) allows AddAllAtoms
       to attempt more জোর (forceful) methods if simple ones fail.
     */
    fprintf (stderr, "Add atoms, pass number 1\n");
    if (i = AddAllAtoms (newnAtoms,newatoms,atomtypes, FALSE)) { /* Initial pass */
	while (TRUE) {
	    fprintf (stderr, "Add atoms, another pass\n");
 	    j = AddAllAtoms (newnAtoms,newatoms,atomtypes, FALSE); /* Standard pass */
 	    k = AddAllAtoms (newnAtoms,newatoms,atomtypes, TRUE);  /* Final (forceful) pass */
	    if (k == 0 || k == i) break; /* If no atoms placed (k==0) or no change from previous full cycle (k==i) */
	    i = k; /* Store the number of unplaced atoms after the forceful pass for next cycle comparison */
	}
	if ( k>0 ) { /* If there are still atoms without coordinates after iterations */
	    fprintf (stderr, "REFORMAT: **** ERROR **** Not all atoms have coordinates\n");
        error = ERR_ATOM_COORDS;
        /* No exit here, proceed to free memory and then return error */
	}
    }

    if (error == 0) { /* Only read params if no errors so far */
        if (readParam (fileparms, newnAtoms, newatoms) != 0) {
            /* Error message already printed by readParam */
            error = ERR_GENERAL; /* Or a more specific error from readParam */
        }
    }

    /* Print atoms and coordinates as roughly PDB */
    (void) WriteXYZ (outfile, newnAtoms, newatoms, residues);

    /* Free allocated memory */
    free(residues);
    free(newatoms);
    free(atomtypes);
    free(restypes);

    if (outfile) fclose(outfile); /* Close the output file if it was opened */
    return error; /* Return 0 on success, or the error code */
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  function ForwardOneBond()
 *    Find the Forwardchain from the info in "atoms[]" and in "atom types[]"
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int ForwardOneBond(int from, PDB *atoms, ATOMTYPE *atomtypes, int nAtoms)
{
int j;
/* int found=FALSE; */ /* Unused */
char *forwardname;
ATOMTYPE *anatomtype;

    anatomtype = atomtypes + atoms[from].key_dict;
    forwardname = anatomtype->forward;
	if (EQUAL(forwardname,"NOT")) return (-1);
    for (j=from+1; j<nAtoms; j++) {
	if (EQUAL(forwardname,atoms[j].atomName)) return (j);
    }
    fprintf(stdout,"REMARK ERROR - FORWARDCHAIN: none found with name %s starting from atom %d\n",
	forwardname,from);
    return (-1);
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  function BackOneBond()
 *    Find the backchain from the info in "atoms[]" and in "atom types[]"
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int BackOneBond(int from, PDB *atoms, ATOMTYPE *atomtypes)
{
int j;
/* int found=FALSE; */ /* Unused */
char *backname;
ATOMTYPE *anatomtype;

    anatomtype = atomtypes + atoms[from].key_dict;
    backname = anatomtype->back;
    if (EQUAL(backname,"NOT")) return (-1);
    for (j=from-1; j>=0; j--) {
	if (EQUAL(backname,atoms[j].atomName)) return (j);
    }
    fprintf(stdout,"REMARK ERROR - BACKCHAIN: none found with name %s starting from atom %d\n",
	backname,from);
    return (-1);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  function readDict
 *    input a file with atom types = a dictionary of atoms in residues
 *    Returns 0 on success, error code on failure.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int readDict (char *filename, int *numat, int *numrestyp, ATOMTYPE **atomtypes, RESTYPE **restypes)
{
FILE *in_file;
ATOMTYPE *anatomtype;
RESTYPE *arestype;
char recdName[10],resname[5];
char term[5],nterm[5],cterm[5];
int iat,i,j;
char line[256];
char *what;

if (*filename == NULLCHAR) {
    what = getenv ("DOWSER");
    if (what) {
        strncpy(filename, what, 256 - 1); /* Assuming filename points to a buffer of 256 */
        filename[256-1] = NULLCHAR;
    }
    else {
	fprintf (stderr,"REMARK ERROR: must first set environment variable 'DOWSER'\n");
	return ERR_DOWSER_ENV_DICT;
    }
    strncat(filename, "/DATA/atomdict.db", 256 - strlen(filename) - 1); /* Assuming filename points to a buffer of 256 */
}
if (!(in_file = fopen (filename,"r"))) {
    fprintf (stderr,"REMARK ERROR: cannot open file %s\n",filename);
    return ERR_OPEN_ATOMDICT;
}

*atomtypes = (ATOMTYPE *) malloc (MAXTYPES * sizeof(ATOMTYPE));
*restypes = (RESTYPE *) malloc (MAXTYPES * sizeof(RESTYPE));

/* The dictionary file (atomdict.db) typically contains:
   - "RESIDUE" lines defining residue names and their N-terminal/C-terminal variants (e.g., "RESIDUE ALA TERM NH3 COO").
   - "ATOM" lines defining atoms within each residue, their connectivity (back, forward),
     ideal geometry (bond, angle, dihedral), charge, and type (e.g., "ATOM ALA  N    C    CA    1.320  114.0  180.0  -0.280    N").
*/

*numat=0; *numrestyp=0;
anatomtype = *atomtypes;
arestype = *restypes;
iat=0;
while (1) {
    fgets (line, 100, in_file); if (feof(in_file)) break;
    sscanf (line, "%s %s", recdName, resname);
    if (EQUAL(recdName,"residue")) {
#ifdef TEST
	fprintf(stderr,"%s\n",resname);
#endif
	term[0]=NULLCHAR;
	sscanf (line, "%s %s %s %s %s", recdName, resname, term, nterm, cterm);
	strncpy(arestype->resname,resname, sizeof(arestype->resname) - 1);
    arestype->resname[sizeof(arestype->resname)-1] = NULLCHAR;
	arestype->begin=iat;
	if (EQUAL(term,"TERM")) {
	    strncpy(arestype->nterminus,nterm, sizeof(arestype->nterminus) - 1);
        arestype->nterminus[sizeof(arestype->nterminus)-1] = NULLCHAR;
	    strncpy(arestype->cterminus,cterm, sizeof(arestype->cterminus) - 1);
        arestype->cterminus[sizeof(arestype->cterminus)-1] = NULLCHAR;
	}
	else *(arestype->nterminus) = *(arestype->cterminus) = NULLCHAR;
    if (*numrestyp >= MAXTYPES) {
        fprintf(stderr, "REMARK ERROR: Number of residue types exceeds MAXTYPES (%d) in atomdict file.\n", MAXTYPES);
        fclose(in_file);
        return ERR_MAXTYPES_EXCEEDED;
    }
	(*numrestyp)++; arestype++;
    }
    if (EQUAL(recdName,"atom")) {
    if (*numat >= MAXTYPES) {
        fprintf(stderr, "REMARK ERROR: Number of atom types exceeds MAXTYPES (%d) in atomdict file.\n", MAXTYPES);
        fclose(in_file);
        return ERR_MAXTYPES_EXCEEDED;
    }
	sscanf (line, "%s %s %s %s %s %f %f %f %f %s",
	    recdName,
	    anatomtype->resname, anatomtype->name, anatomtype->back, anatomtype->forward,
	    &anatomtype->bond, &anatomtype->angle, &anatomtype->dihedral,
	    &anatomtype->charge, anatomtype->type);
	anatomtype->index = iat;
	(arestype-1)->end=iat;
	(*numat)++; anatomtype++; iat++;
    }
}
arestype->end=iat;

#ifdef TEST
arestype = *restypes;
for (i=0; i<*numrestyp; i++) {
    fprintf(stdout,"RES: %s %d %d\n",arestype->resname,arestype->begin,arestype->end);
    for (j=arestype->begin; j<=arestype->end; j++) {
	anatomtype = *atomtypes+j;
	fprintf(stdout,"ATOM: %s %s\n", anatomtype->name,anatomtype->type);
    }
    arestype++;
}
#endif
fclose(in_file);
return 0; /* Success */
} /* end readDict() */

/*
ATOM    349  OD1 ASN    43      29.384  10.062  -2.624  1.00  4.58      1BPI 485
******+++++ ****?*** +****?   ********++++++++********OOOOOOBBBBBB ???  ****
12345678901234567890123456789012345678901234567890123456789012345678901234567890
*/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  WriteXYZ() : write the new pdb file with LJ parameters and charges
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void WriteXYZ(FILE *outFileP, int num, PDB *atoms, RESIDUE *residues)
{
int i, j, rnum;
PDB *anatom;
char aname[5],achar;
char *what;
int backchain;

    for (i=0;i<num;i++) {
	anatom = atoms +i;
	if (!anatom->ftNote) continue;
	rnum = anatom->newResNo;
	strncpy(aname,anatom->atomName, sizeof(aname) - 1);
    aname[sizeof(aname)-1] = NULLCHAR;
	if (strlen(aname) == 4) { /* This check might be problematic if aname was truncated */
	    achar = aname[3];
	    aname[3] = NULLCHAR;
	}
	else achar = ' ';
	what=atomtypes[anatom->key_dict].back;
	if (EQUAL(what,"NOT")) backchain=0;
	else {
	    for (j=i-1;j>=0;j--) if (EQUAL(atoms[j].atomName,what)) break;
	    if (j==-1) {
		backchain=0;
		fprintf(stderr,"Bad backchain atom %d %s, %s\n",i,anatom->atomName,what);
	    }
	    else backchain = i-j;
	}
	fprintf(outFileP, 
	    "ATOM  %5d %c%-4s%-3s %5d    %8.3f%8.3f%8.3f  %4.2f%6.2f %7.3f %7.2f %7.1f %6.2f %-4s %3d\n",
	    i+1, achar, aname, residues[rnum].name,
	    rnum + 1, anatom->XX, anatom->YY, anatom->ZZ,
	    anatom->occupancy, anatom->tempFactor,
	    anatom->LJ_c, anatom->LJ_a, anatom->LJ_b, anatom->ms_rad,
	    atomtypes[anatom->key_dict].type,backchain);
    }
    fprintf (stdout,"* Reformatted and extended pdb file contains %d atoms\n", num);
}
/*
ATOM      2  CA  ACE     1      -2.971  -0.162   1.364
123456789012345678901234567890123456789012345678901234567890
*/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  readParam() : 
 *     read the LJ parameters and assign them to atoms
 *    Returns 0 on success, error code on failure.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int readParam(char *filename, int nAtoms, PDB *atoms)
{
PDB *anatom;
char line[256];
char *what;
FILE *ffl;
/* int eof = FALSE; */ /* Unused */
char type[10];
REAL a, b, ms_rad;
int i;

    if (*filename == NULLCHAR) {
	what = getenv ("DOWSER");
	if (what) {
        strncpy(filename, what, 256 - 1); /* Assuming filename points to a buffer of 256 */
        filename[256-1] = NULLCHAR;
    }
	else {
	    fprintf (stderr,"REMARK ERROR: must first set environment variable 'DOWSER'\n");
	    return ERR_DOWSER_ENV_PARAM;
	}
	strncat(filename, "/DATA/atomparms.db", 256 - strlen(filename) - 1); /* Assuming filename points to a buffer of 256 */
    }

    if (!(ffl = fopen (filename,"r"))) {
	fprintf (stderr,"REMARK ERROR: cannot open file %s\n",filename);
	return ERR_OPEN_ATOMPARMS;
    }

    while (TRUE) {
	fgets(line, 100, ffl);
	if (feof(ffl)) break;
	if (!strncmp(line,"TYPE",4)) {
	    sscanf (line+7, "%s", type);
	    sscanf (line+10, "%f", &a);
	    sscanf (line+19, "%f", &b);
	    sscanf (line+29, "%f", &ms_rad);
	    anatom = atoms;
	    for (i=0; i<nAtoms; i++) {
		if (EQU(anatom->atomType,type)) {
		    anatom->LJ_a = a;
		    anatom->LJ_b = b;
		    anatom->ms_rad = ms_rad;
		    anatom->ftNote = 2;
		}
		anatom ++;
	    }
	}
    }
    fclose (ffl);
    return 0; /* Success */
}
/*
REMARK atomtype LJ-a LJ-b    MS-radius
1234567890123456789012345678901234
TYPE   N      49.36   1300.0  2.40
TYPE   N      49.36   1300.0
012345678901234567890123456789
*/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *   FindSSBonds ( )
 *     Use a distance criterion 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void FindSSBonds (FILE *outFileP, PDB *atoms, RESIDUE *residues, int nRes, RESTYPE *restypes, int nRestype)
{
int i,i1,i2,j,k;
REAL xs1[3],xs2[3];
REAL ds, ds2;
static int css_type=-1;
/* char *what; */ /* Unused */

if (css_type < 0) {
    for (css_type=0; css_type < nRestype; css_type++ )
	if (EQUAL("CSS",restypes[css_type].resname)) break;
    if (css_type==nRestype) {
	fprintf(stderr,"REMARK ERROR: Residue type CSS not found in the type dictionary\n");
	fprintf(stderr,"REMARK ERROR: Will not be able to handle SS-bonds\n");
    }
}
if (css_type==nRestype) return;

  /* locate disulfide bridges and change residue names */
  for (i1=0;i1<nRes;i1++) {
      if (j=CysSG (atoms,i1,residues,xs1)) {
	  /* look for another CYS residue */
	  for (i2=i1+1; i2<nRes; i2++) {
	      if (k=CysSG (atoms,i2,residues,xs2)) {
		  ds2=0.;
		  for (i=0;i<3;i++) { ds = xs2[i] - xs1[i]; ds2 += ds*ds; }
		  if ( ds2 <= SSBONDSQ ) {
		      /* strcpy(residues[i1].name,"CSS"); */
		      residues[i1].type=css_type;
		      /* strcpy(residues[i2].name,"CSS"); */
		      residues[i2].type=css_type;
fprintf(outFileP,"REMARK SS-BOND residues %d and %d, ds2 = %f\n",i1+1,i2+1,ds2);
		      break;
		  }
	      }
	  }
      }
  }
} /* end FindSSBonds() */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  int CysSG (atoms,resnum,atoms,xs1)
 *  return atom index of SG if residue is a CYS, else 0
 *  and the coordnates of this atom
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int CysSG (PDB *atoms, int resnum, RESIDUE *residues, REAL *xs1)
{
int j=0;
  if (EQUAL(residues[resnum].name,"CYS")) {
      if (j=LocateInResidue(atoms,resnum,residues,"SG")) {
	  xs1[0] = atoms[j].XX;
	  xs1[1] = atoms[j].YY;
	  xs1[2] = atoms[j].ZZ;
      }
  }
  return (j);
} /* end CysSG() */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  int LocateInResidue(atoms,resnum,residues,name)
 *     return index of atom with given name
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int LocateInResidue(PDB *atoms, int resnum, RESIDUE *residues, char *name)
{
int j;
    for (j=residues[resnum].begin0;j<=residues[resnum].end0;j++) {
      if (EQUAL(atoms[j].atomName,name)) return (j);
    }
    return (0);
} /* end LocateInResidue() */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  FindChainBreaks(atoms,residues,nRes,restypes)   
 *     Use a distance criterion 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void FindChainBreaks(FILE *outFileP, PDB *atoms, RESIDUE *residues, int nRes, RESTYPE *restypes, ATOMTYPE *atomtypes)
{
int ires,i;
int j1,j2;
REAL dlink, dlink2;
int num_breaks=0;
RESTYPE *arestype;
char *aname,*bname;
int type,begin;

/* chain termination for first and last molecule */
arestype = restypes + residues[0].type;
strncpy(residues[0].terminus,arestype->nterminus, sizeof(residues[0].terminus) - 1);
residues[0].terminus[sizeof(residues[0].terminus)-1] = NULLCHAR;
residues[0].chainend=1;
arestype = restypes + residues[nRes-1].type;
strncpy(residues[nRes-1].terminus,arestype->cterminus, sizeof(residues[nRes-1].terminus) - 1);
residues[nRes-1].terminus[sizeof(residues[nRes-1].terminus)-1] = NULLCHAR;
residues[nRes-1].chainend=2;

for (ires=1;ires<nRes;ires++) {
    /* first atom in the residue according to the dictionary */
    type = residues[ires].type;
    begin = restypes[type].begin;
    aname = atomtypes[begin].name;

    j2 = LocateInResidue(atoms,ires,residues,aname);
    bname = atomtypes[begin].back;
    if (EQUAL(bname,"NOT")) j1=0; /* always a break */
    else j1 = LocateInResidue(atoms,ires-1,residues,"C");

    /* check distance of link */
    dlink2=0;
    if (j1 && j2) {
	dlink = atoms[j1].XX - atoms[j2].XX;
	dlink2 += dlink*dlink;
	dlink = atoms[j1].YY - atoms[j2].YY;
	dlink2 += dlink*dlink;
	dlink = atoms[j1].ZZ - atoms[j2].ZZ;
	dlink2 += dlink*dlink;
	if (dlink2 <= NCBONDSQ) continue; /* it's a bond */
    }

    fprintf (outFileP,"REMARK CHAIN BREAK at residue %d, C-N 'bondlength' = %8.2f\n", ires+1, sqrt(dlink2));
    num_breaks++;

    arestype = restypes + residues[ires].type;
    residues[ires].chainend=1;
    strncpy(residues[ires].terminus,arestype->nterminus, sizeof(residues[ires].terminus) - 1);
    residues[ires].terminus[sizeof(residues[ires].terminus)-1] = NULLCHAR;
    arestype = restypes + residues[ires-1].type;
    residues[ires-1].chainend=2;
    strncpy(residues[ires-1].terminus,arestype->cterminus, sizeof(residues[ires-1].terminus) - 1);
    strncpy(residues[ires-1].terminus,arestype->cterminus, sizeof(residues[ires-1].terminus) - 1);
    residues[ires-1].terminus[sizeof(residues[ires-1].terminus)-1] = NULLCHAR;
}

if (!num_breaks) fprintf(outFileP,"REMARK NO CHAIN BREAKS\n");

} /* end FindChainBreaks() */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *    void MergeTerminus
 *        (countonly,atoms,aresidue,nRestype,restypes,atomtypes,
 *         mergeatomtype,mergerestype)
 * Returns 0 on success, error code on failure.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int MergeTerminus
    (FILE *outFileP, int countonly, PDB *atoms, RESIDUE *aresidue, int nRestype, RESTYPE  *restypes, ATOMTYPE *atomtypes, ATOMTYPE *mergeatomtype, RESTYPE *mergerestype)
{
char *what; /* Temporary string pointer, often to atom names */
RESTYPE *arestype, *brestype, *res_restype, *term_restype;
int found, i, j;
int keytype;
ATOMTYPE *anatomtype, *batomtype;

    /* find the terminus residue type (e.g., "NH3", "COO") and the original residue type (e.g., "ALA")
       from the global restypes array. */
    what=aresidue->terminus; /* Name of the terminus type, e.g., "NH3" */
    for (i=0;i<nRestype;i++) {
	arestype = restypes+i;
	if (EQUAL(arestype->resname,aresidue->terminus)) term_restype=arestype; /* Found the RESTYPE for the terminus */
	if (EQUAL(arestype->resname,aresidue->name))     res_restype=arestype;  /* Found the RESTYPE for the original residue name */
    }

    if (countonly) { /* First pass: just count extra atoms needed for the terminus */
        /* Iterate through atoms in the terminus dictionary type (e.g., atoms of NH3) */
	for (i=term_restype->begin;i<=term_restype->end;i++) {
	    what=atomtypes[i].name; /* Name of an atom in the terminus type */
	    found=FALSE;
            /* Check if this terminus atom is already defined in the original residue type */
	    for (j=res_restype->begin;j<=res_restype->end;j++) {
		if (EQUAL(atomtypes[j].name,what)) {
		    found = TRUE; break; /* Atom already exists in original residue type */
		}
	    }
	    if (!found) aresidue->numextra++; /* If not found, it's an extra atom to be added */
	}
	if (aresidue->numextra)
	fprintf (outFileP,"REMARK Need %d extra atoms in residue %s %d\n",
	    aresidue->numextra, aresidue->name, aresidue->index);
    }
    /* Second pass (actual merge) or if not countonly:
       Construct the list of atom types for the modified (terminal) residue.
       The mergeatomtype array will hold the final list of ATOMTYPEs for this residue.
       mergerestype will describe this newly defined merged residue type.
    */
    else {
	keytype=0; /* Index for mergeatomtype array */

	/* If it's an N-terminus (aresidue->chainend==1):
	   First, add all atoms from the N-terminus definition (e.g., NH3 atoms). */
	if (aresidue->chainend==1) { /* insert N-terminus */
	    for (i=term_restype->begin;i<=term_restype->end;i++) {
        if (keytype >= MAXMERGE) {
            fprintf(stderr, "REMARK ERROR: Merged atom count exceeds MAXMERGE (%d) in residue %s %d for N-terminus.\n", MAXMERGE, aresidue->name, aresidue->index);
            return ERR_MAXMERGE_EXCEEDED;
        }
		mergeatomtype[keytype++] = atomtypes[i]; /* Add N-terminus atom type */
    }
	} /* end chain begin */

	/* Next, add atoms from the original residue type (e.g., ALA),
	   but only if an atom with the same name wasn't already added from the N-terminus. */
	arestype = restypes + aresidue->type; /* Original residue type from PDB (e.g. ALA, CYS) */
	for (i= arestype->begin; i<=arestype->end; i++) {
	    anatomtype = atomtypes + i; /* Current atom from original residue type */
	    what = anatomtype->name; 
            /* Check if this atom name already exists in mergeatomtype (from N-terminus) */
	    for (batomtype=mergeatomtype;batomtype<mergeatomtype+keytype;batomtype++)
		if (EQUAL(batomtype->name,what)) goto SKIP_ME; /* Atom already added, skip */
        if (keytype >= MAXMERGE) {
            fprintf(stderr, "REMARK ERROR: Merged atom count exceeds MAXMERGE (%d) in residue %s %d.\n", MAXMERGE, aresidue->name, aresidue->index);
            return ERR_MAXMERGE_EXCEEDED;
        }
	    mergeatomtype[keytype++] = *anatomtype; /* Add original residue atom type */
	    SKIP_ME: ;
	} /* end residue */

	/* If it's a C-terminus (aresidue->chainend==2):
	   Add atoms from the C-terminus definition (e.g., COO atoms).
	   If an atom with the same name already exists (e.g. from original residue, or N-term if it's a single AA),
	   the C-terminus definition *overwrites* the previous one. This is important for atoms like 'O'.
	*/
	if (aresidue->chainend==2) {
	    for (i=term_restype->begin;i<=term_restype->end;i++) {
		anatomtype = atomtypes + i; /* Current atom from C-terminus type */
		what = anatomtype->name; 
                /* Check if this atom name already exists in mergeatomtype */
		for (batomtype=mergeatomtype;batomtype<mergeatomtype+keytype;batomtype++) {
		    if (EQUAL(batomtype->name,what)) {
			*batomtype = *anatomtype; /* Atom exists, overwrite its specs with C-terminus version */
			goto SKIP_ME2; /* Go to next C-terminus atom */
		    }
		}
        /* If atom was not found in mergeatomtype, add it as a new entry */
        if (keytype >= MAXMERGE) {
            fprintf(stderr, "REMARK ERROR: Merged atom count exceeds MAXMERGE (%d) in residue %s %d for C-terminus.\n", MAXMERGE, aresidue->name, aresidue->index);
            return ERR_MAXMERGE_EXCEEDED;
        }
		mergeatomtype[keytype++] = *anatomtype; /* Add C-terminus atom type */
		SKIP_ME2: ;
	    }
	} /* end chainend */
    }
    strncpy(mergerestype->resname,aresidue->name, sizeof(mergerestype->resname) - 1);
    mergerestype->resname[sizeof(mergerestype->resname)-1] = NULLCHAR;
    mergerestype->begin=0;
    mergerestype->end= keytype-1;
    mergerestype->numat= keytype;

    return 0; /* Success */

} /* end MergeTerminus */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *    int ScanDictVsPDB(nRes,residues,nRestype,restypes)
 *    compare the sequence in the PDB wioth the available residue types
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int ScanDictVsPDB(int nRes, RESIDUE *residues, int nRestype, RESTYPE *restypes)
{
int i,j;
int error=FALSE;
char *what;

    for (i=0; i<nRes; i++) {
	what = residues[i].name;
	for (j=0; j < nRestype; j++ ) {
	    if (EQUAL(what,restypes[j].resname)) { residues[i].type = j; break; }
	}
	if (j==nRestype) {
	    fprintf(stderr,"REMARK ERROR: Residue type %s not found in the type dictionary\n",what);
	    error=TRUE;
	}
    }
    return (error);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *    int AddAllAtoms()
 *        Compute the coordinates of the missing atoms 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int AddAllAtoms (int nAtoms, PDB *atoms, ATOMTYPE *atomtypes, int final_pass)
{
int i, fault = 0; /* fault counts atoms whose coordinates could not be determined */
int refatom, backatom, backback;
ATOMTYPE *anatomtype, *bnatomtype, *cnatomtype;

for (i=0;i<nAtoms; i++) {
    if (!atoms[i].ftNote) {
	atoms[i].occupancy = 0.;
	atoms[i].tempFactor = 0.;
	anatomtype = atomtypes + atoms[i].key_dict; /* ATOMTYPE of the current atom 'i' to be placed */

	/* Strategy for placing atom 'i':
	   Ideally, uses its direct bonded neighbors (back, forward) and angles/dihedrals from dictionary.
	   Atoms are defined by internal coordinates: bond length from atom 'c', angle a-c-b, dihedral a-c-b-d.
	   'i' is 'a'. 'backatom' is 'c'. 'backback' is 'd'. 'refatom' (from forward chain) is 'b'.
	*/

	/* Try to find the atom 'c' (atom 'i' is bonded to, going back in sequence) */
	backatom = BackOneBond (i,atoms,atomtypes);

	/* Case 1: Atom 'i' has no back-bonded atom (e.g., N-term) OR it's the final pass.
	   In this 'root atom' or 'final pass' scenario, try to build using forward chain atoms as reference.
	   This means atom 'i' is 'd', 'backatom' (from forward) is 'b', 'backback' (from forward) is 'c', 'refatom' (from forward) is 'a'.
	   This is effectively reversing the reference frame.
	*/
	if (backatom<0 || final_pass) { 
	    backatom = ForwardOneBond (i,atoms,atomtypes,nAtoms); /* This will be 'b' for Add1Atom */
	    if (backatom<0) continue; /* Cannot find reference, skip this atom for now */
	    backback = ForwardOneBond (backatom,atoms,atomtypes,nAtoms); /* This will be 'c' for Add1Atom */
	    if (backback<0) continue;
	    refatom =  ForwardOneBond (backback,atoms,atomtypes,nAtoms); /* This will be 'd' for Add1Atom (planar ref)*/
	    if (refatom<0) continue;

        /* Get ATOMTYPE definitions for the reference atoms to extract geometry parameters */
	    ATOMTYPE *atom_i_type = atomtypes + atoms[i].key_dict; /* for atom to be placed */
	    ATOMTYPE *atom_b_type = atomtypes + atoms[backatom].key_dict; 
        ATOMTYPE *atom_c_type = atomtypes + atoms[backback].key_dict;
        /* For Add1Atom(atoms, d, c, b, a, bond_ca, angle_dcb, dihedral_dcba)
           We are placing 'i' (acting as 'a' in Add1Atom's frame, but it's the 4th atom in chain d-c-b-a)
           Here, 'refatom' is d, 'backback' is c, 'backatom' is b.
           Bond length is c-b (atom_b_type->bond, assuming 'forward' implies 'bond' to previous)
           Angle is d-c-b (atom_c_type->angle)
           Dihedral is d-c-b-i (atom_i_type->dihedral, relative to its forward connection)
           This part of the logic might need careful review of how Add1Atom expects its args
           when building "backwards" or from a different reference.
           The current Add1Atom call uses:
           Add1Atom(atoms, d=backback, c=backatom, b=refatom, a=i, ...)
           This seems to be trying to place 'i' using 'refatom' as the third atom in a chain.
           The geometry parameters are taken from types of refatom, backback, backatom.
        */
        /* The following parameters seem to be for placing 'i' based on 'refatom', 'backback', 'backatom' as d,c,b */
	    ATOMTYPE *type_of_d =  atomtypes + atoms[refatom].key_dict;  /* d */
	    ATOMTYPE *type_of_c =  atomtypes + atoms[backback].key_dict; /* c */
	    ATOMTYPE *type_of_b =  atomtypes + atoms[backatom].key_dict;   /* b */
	    if (Add1Atom (atoms,refatom, backback, backatom, i, /* d, c, b, a (atom to place) */
	      type_of_b->bond, /* bond c-b */
	      (PI/180.)*type_of_c->angle, /* angle d-c-b */
	      (PI/180.)* type_of_d->dihedral) ) { /* dihedral for atom 'd' when connected to c,b,a. This seems off. */
		atoms[i].ftNote = TRUE; /* Mark as coordinates calculated */
	    }
	    continue; /* Move to next atom */
	}

	/* Case 2: Atom 'i' has a valid 'backatom' (atom 'c'). This is the standard case. */
	backback = BackOneBond (backatom,atoms,atomtypes); /* Atom 'd' */
	refatom = ForwardOneBond (backatom,atoms,atomtypes,nAtoms); /* Atom 'b' (forward from 'c') */
	
    /* If 'refatom' ('b') is a valid atom (not 'H', as H might be placed later) and further in sequence than 'i' (avoiding loops) */
	if (refatom > i && atoms[refatom].atomName[0] != 'H') {
	    bnatomtype =  atomtypes + atoms[refatom].key_dict; /* Type of atom 'b' */
        /* Add1Atom(atoms, d, c, b, a, bond_ca, angle_dcb, dihedral_dcba)
           Here, a=i, c=backatom, d=backback, b=refatom.
           bond_ca is anatomtype->bond. angle_dcb is anatomtype->angle.
           dihedral_dcba is (anatomtype->dihedral - bnatomtype->dihedral) - this is a relative dihedral.
        */
	    if (Add1Atom (atoms,backback, backatom, refatom, i,
	      anatomtype->bond, (PI/180.)*anatomtype->angle,
	      (PI/180.)* (anatomtype->dihedral - bnatomtype->dihedral)) )
	      atoms[i].ftNote = TRUE;
	}
	else {
	    /* Case 3: Standard 'backatom' exists, but 'refatom' (forward) is not suitable.
	       Use a chain of three preceding atoms: back(backback) -> backback -> backatom -> i
	       Here, refatom becomes back(backback) ('e'), backback is 'd', backatom is 'c'.
	    */
	    if (backback>=0) refatom = BackOneBond (backback,atoms,atomtypes); /* Atom 'e' */
	    if (backback>=0 && refatom >=0) { /* If d and e exist */
            /* Add1Atom(atoms, d, c, b, a, bond_ca, angle_dcb, dihedral_dcba)
               Here, a=i, c=backatom, d=backback, b=refatom (which is actually 'e').
               This reuses refatom variable for the third point in the reference chain.
               So, effectively: Add1Atom(atoms, backback, backatom, refatom_as_e, i, ...)
            */
	      if (Add1Atom (atoms,backback, backatom, refatom, i,
		anatomtype->bond, (PI/180.)*anatomtype->angle, (PI/180.)*anatomtype->dihedral) )
		atoms[i].ftNote = TRUE;
	    }
	}
    }
}

for (i=0;i<nAtoms; i++) {
    fault += ( 1 - atoms[i].ftNote );
}
fprintf (stderr, "%d atoms not found\n", fault);
return fault;
} /* end AddAllAtoms() */
