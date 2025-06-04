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
 *  AmideFlip.c
 *  by Philip M-H Kim, modified by Dave Cavanaugh and Jan Hermans
 *  Resolve amide ambiguities by flipping them to lowest energy configuration.  Also,
 *  handle sets of connected amides (e.g. GLN - ASN pairs).
 *
 *  Input:
 *    argv[1] - input pdb structure (should be output from reformatPDB)
 *  Output:
 *    argv[2] - output pdb structure
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

#define SHOWPROGRESS
#define SHOWGROUPS
#undef DEBUG

#define NB_CUTOFF_SQ 1.e6  /* if 2 atoms are farther apart, the NB energy is not calculated */
#define GROUP_CUTOFF  7    /* if 2 amides are within this distance, they are grouped */

#define AMIDE(i) pro[i].key_dict
#define MOVE(i)  pro[i].newResNo
#define NULL_CHAR '\0'

double sum_e_el,sum_e_lj;

int amNum;               /* Number of total amides */ 
int connum;              /* Number of connectivities in current amidegroup */
PDB *pro;                /* Protein Data */

typedef struct { int iat,jat; } Pair;
Pair *NBPairs; 		/* pair-list */
int NumberNBPairs; 

int connect[100];        /* Amidenumbers that are connected to current amide */
int **resorder;          /* Amides that are connected to each other are written in one row */
int **atoms_in_am;         /* Atom numbers for diefferent amide atoms */
int **distance;          /* Distancematrix, 1 when amides are connected, 0 otherwise */
int numAtms;       /* Number of atoms, i,j indices */

int resgroup_num;      /* Number of amidegroups */

int SetAmideFlag(PDB *pro_param);
int GroupAmides(void);
void FlipOneAmide(int flipme);
REAL AmideEnergy(void);
int AmidePairList(int *resorder_i);
void AtomPairs(int central_atom, int *numpairs);
void Search(int residnum, int *prior);
int check(int number);
int testcon(int rescon);
int extbit(int string, int bit);

extern void readPDB10(char *filename, int *num_atoms, PDB **atoms_ptr);
extern void initVDW(char *filename, int num_atoms, PDB *atoms);
extern void WriteXYZ(FILE *outfile, int num, PDB *atoms);
extern REAL DistSq(REAL r1[3], REAL r2[3]);

/*  Main */
int main(int argc, char *argv[])
{ 
  
/* Give listing of what variables represent */
  char flag_for_fl[10];
  int amNum_group;       /* Numbers of amides in current amidegroup (resgroup) */
  int res1, res2;        /* Residuenumbers for distancematrix */
  int resnum;

  int k, configNum;      /* k: index, configNum: Number of different configurations in current amidegroup */
  int second_best;       /* For calculation of energy difference: The configuration of amidegroup in which 
			    the current (flipped) amide was unflipped ! */
  int am_pdb_num;        /* Number of amideresidue in PDB-file (for reading purposes) */
  int i,j;      	/* Number of atoms, i,j indices */
  int num;               /* The Number of the current amide in its group */
  REAL En_for_group;     /* Energydifference for whole group */

  REAL *Energy_flip;    /* energies for flipped alternative arrangements of a group of residues */
  REAL min_energy;	/* min. energy value for current group */
  int *min_en_conf;	/* min. energy configuration for all groups */

  int count;             /* count atom numbers */
  int amide1,amide2;     /* Amide atom numbers */
  int atoms_in_res;      /* Number of atoms in current residue */
  int max_am_group;      /* Maximum Number of amides in group */
  REAL atomdis;          /* Distance between amides (for Distancematrix calculation) */

  int max_group_size;	/* max. number of amides in a connected group */
  int have_flipped;

  FILE *fp;

  numAtms = 0; /* Initialize variables to 0 */
  num = resgroup_num = 0;
  count =0;
  amide1 = amide2 = 0;
  atoms_in_res = 0;
  amNum_group =0;
  resnum = connum = k =configNum =0;
  NumberNBPairs =0;

  for (i=0; i<100; i++) connect[i]=0;
  

  if (argc != 3) {
      fprintf(stderr,"USAGE: AmideFlip input.pdb output.pdb\n");
      exit(1);
  }
  readPDB10(argv[1],&numAtms,&pro);

/* read nonbonded parameters */
  initVDW (argv[1], numAtms, pro);
  
  fprintf(stdout,"* FLIP: Number of atoms in Protein: %d \n",numAtms);

  /* * * * * * * Assign the AMIDE(i) flags * * * * * */
  /* if not amide, flag=0; if one of C-C(=O)-NH2 flag=index of this amide; 
    if preceding atom, flag = -index of the amide */
    amNum = SetAmideFlag(pro);
    fprintf(stdout,"* FLIP: Number of Amides = %d\n",amNum); 

  /* group the amides by proximity */
  amNum_group = GroupAmides();
  min_en_conf = (int *) malloc (amNum_group*sizeof(int));

  max_group_size=0;
  for (i=0; i<resgroup_num;i++) {
    j=0;
    while (resorder[i][j] != -1) {
	j++;
	if (j>max_group_size) max_group_size=j;
	if (j>31) {
	    fprintf(stderr,"WARNING, more than 32 amides in one group, can't handle that !");
	    exit(1);
	}
    }
  } 

  configNum=pow(2,max_group_size);
  Energy_flip = (REAL *) malloc (configNum * sizeof (REAL));
  if (!Energy_flip) {
      fprintf (stderr,"** Unable to allocate memory for Energy_flip\n");
      exit (1);
  }

/*************************************************************************************************/
/* The energies can now be evaluated */

    /* consider each group of side chains in turn */
    for (i=0; i<resgroup_num; i++) {
	min_en_conf[i]=0;
	/* how many residues in this group */
	amNum_group=0;
	while (resorder[i][amNum_group]!=-1)  amNum_group++;
	if (amNum_group == 0) continue;

	configNum=pow(2,amNum_group); /* number of flip configurations */
#ifdef SHOWPROGRESS
	fprintf (stderr,"* Config.  Total, el., lj energies: for group number %d\n",i);
#endif
	/* pairlist for this group */
	NumberNBPairs = AmidePairList(resorder[i]);

	for (j=0; j<configNum; j++) {
	    /* flip the amides according to "1" bits in k */
	    for (k=0; k<amNum_group; k++) {
		if (extbit(j,k)==1) { FlipOneAmide(resorder[i][k]); }
	    } /* all flipped that need to be flipped */

#ifdef SHOWPROGRESS
	    fprintf(stderr,"*   %3d,",j);
#endif
	    /* calculate the energy for for this configuration */
	    Energy_flip[j] = AmideEnergy();

	    /* must restore "all 0 configuration" */
	    for (k=0; k<amNum_group; k++) {
		if (extbit(j,k)==1) FlipOneAmide(resorder[i][k]);
	    }
	} /* end j-loop over possible flips of one group */

	min_energy=1.e10; 
	for (j=0;j<configNum;j++) {
	    if (Energy_flip[j] < min_energy) {
		min_energy=Energy_flip[j]; min_en_conf[i]=j;
	    }
	}
#ifdef DEBUG
	fprintf(stdout,"* Minimum energy of group %d= %f, in config. %d\n",
	    i,min_energy,min_en_conf[i]);
#endif
      } /* end i-loop over groups of linked side chains */

/* Flip residues to minimum energy configuration: */
  have_flipped = FALSE;
  for (i=0; i<resgroup_num; i++) {
    if ((j=min_en_conf[i]) == 0) continue;
    amNum_group=0;
    while (resorder[i][amNum_group]!=-1) amNum_group++;
    if (!have_flipped) {
	fprintf(stdout,"* Some amides will be flipped to obtain a lower-energy conformation\n");
	have_flipped = TRUE;
    }
    printf("*     Flip these amides for group %d:",i);
    for (k=0; k<amNum_group; k++) {
	if (extbit(j,k)==1) {
	    FlipOneAmide((res1=resorder[i][k]));
	    amide1 = atoms_in_am[res1][0];
	    fprintf(stdout," [(%d) %d %s]",res1,pro[amide1].resSeq,pro[amide1].resName);
	}
    }
    fprintf (stdout,"\n");
  }

/* Output resulting structure */
  if (have_flipped) {
    write_out:
    fprintf(stdout,"* FLIP: write modified structure on file %s\n",argv[2]);
    fp = fopen(argv[2],"w");
    if (!fp)  {
	fprintf(stderr,"*** Unable to open %s for output - use stdout instead!\n",argv[2]);
	fp = stdout;
    }
    WriteXYZ (fp,numAtms,pro);
  }
  else fprintf (stdout,"* FLIP: no residues have been flipped, no new structure output\n");
  return 0;
}
/* * * * *  end main * * * * * */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * void FlipOneAmide()
 *    Amide flipping function
 *
 * This performs the folowing calculation:
 *    Project the CO-vector on the CC axis, calculate normal from O-xyz to CC
 *    add 2 * the normal to O-xyz  --> has been flipped ! 
 *    xyz-O += 2 * [ ( <Cam-C,O-Cam> * (Cam-C)/ |Cam-C|^2 ) - (O - Cam)];
 *    Same thing with N, H1 and H2 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void FlipOneAmide(int flipme)
{
    int i,j;
    int i1,j1;
    REAL xyzNew[18];
    REAL vspindle[3];
    REAL vector[18];
    REAL comp;
    REAL projection[6];
    REAL spindle_sq;
    int *atoms;

#ifdef DEBUGX
    fprintf (stdout,"FLIP amide with index %d\n",flipme);
#endif

    atoms = atoms_in_am[flipme];

    /* First assign Corrdinates to temporary Variable */
    i=0;
    for (i1=0; i1<6; i1++) {
	for (j=0;j<3;j++) 
	    xyzNew[i++]=pro[atoms[i1]].xyz[j];
    }


    /* Calculate VECTOR Cam-C */
    spindle_sq=0.;
    for (i1=0; i1<3; i1++) {
	comp=vspindle[i1]= xyzNew[i1+3]-xyzNew[i1];
	spindle_sq += comp*comp;
    }
		 
    /*Calculate VECTORs O-Cam, N-Cam etc.
    Counting from 2 to 6 !!!! b/c 2nd atom is Cam !!*/
    for (i1=2; i1<6; i1++){
	for (j1=0; j1<3; j1++) vector[j1+3*i1]=xyzNew[j1+3*i1]-xyzNew[j1+3];
    }

    for (i1=2; i1<6; i1++) projection[i1]=0.;
    /* Calculate Projection <Cam-C,O/H/N-C)> */
    for (i1=2; i1<6; i1++) {
	for (j1=0; j1<3; j1++) {
	    comp=vspindle[j1] * vector[j1+3*i1];
	    projection[i1] +=comp;
	}
    }

#undef DEBUGFLIP
#ifdef DEBUGFLIP
    fprintf (stderr,"Xyz\n");
    for (i=0; i<18; i++) {
	fprintf (stderr," %8.2f",xyzNew[i]);
	if (i%3 == 2) fprintf (stderr,"\n");
    }
    fprintf (stderr,"\n");
    fprintf (stderr,"Spindle\n");
    for (i=0; i<3; i++) fprintf (stderr," %8.2f",vspindle[i]);
    fprintf (stderr," %8.2f",spindle_sq);
    fprintf (stderr,"\n");
    fprintf (stderr,"Vectors\n");
    for (i=6; i<18; i++) {
	fprintf (stderr," %8.2f",vector[i]);
	if (i%3 == 2) fprintf (stderr,"\n");
    }
    fprintf (stderr,"\n");
    fprintf (stderr,"Projection\n");
    for (i=0; i<6; i++) fprintf (stderr," %8.2f",projection[i]);
    fprintf (stderr,"\n");
#endif

    /* Finally, DO the flipping */

#ifdef DEBUGFLIP
    fprintf(stderr,"xyz and increments\n");
    for (i1=2; i1<6; i1++) {
	for (j1=0; j1<3; j1++) {
	    fprintf(stderr," %8.3f %8.3f,",
	    xyzNew[j1+3*i1],2*(projection[i1]*vspindle[j1]/spindle_sq)-2*vector[j1+3*i1]);
	    xyzNew[j1+3*i1] +=
		2*( projection[i1] * vspindle[j1] / spindle_sq - vector[j1+3*i1] ); 
	} fprintf (stderr,"\n");
    }
#else
    for (i1=2; i1<6; i1++) {
	for (j1=0; j1<3; j1++)
	    xyzNew[j1+3*i1] +=
		2 * ( projection[i1] * vspindle[j1] / spindle_sq - vector[j1+3*i1] ); 
    }
#endif

    /* Copy new coordinates to pro */
    i=6;
    for (i1=2; i1<6; i1++) {
    for (j=0;j<3;j++) 
	pro[atoms[i1]].xyz[j]=xyzNew[i++];
    }

    return;
}

/********************************************************************
 *
 * AmideEnergy() - Calculate total potential energy for the five 
 *                 amide atoms.
 *
 *  returns total energy in kcal/mol
 ********************************************************************/

REAL AmideEnergy(void)
{
int i;
int iat,jat;
REAL c;
REAL a1,b1,c1;          /* parameters for bonded atom */
REAL a2,b2,c2;          /* parameters for amide atom */
double r,r2;            /* distance between protein and water atoms */
REAL ener;              /* electrostatic energy */
REAL e_el;
double r6,r12,e_lj6,e_lj12;
REAL sum_e_el = sum_e_lj = 0.;

    for (i=0; i<NumberNBPairs; i++) {
	iat =NBPairs[i].iat;
	jat =NBPairs[i].jat;

	a1 = pro[jat].LJ_a; b1 = pro[jat].LJ_b; c1 = pro[jat].LJ_c;
	a2 = pro[iat].LJ_a; b2 = pro[iat].LJ_b; c2 = pro[iat].LJ_c;

	/*
	for (j=0;j<3;j++) dist[j] = pro[iat].xyz[j] - pro[jat].xyz[j];
	r2 = 1/(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2]);
	*/

	r2 = 1./DistSq(pro[iat].xyz,pro[jat].xyz);

	if ((c1 == 0) || (c2 == 0))  e_el = 0.;
	else  {
#undef DIST_DEP_DIEL
#ifdef DIST_DEP_DIEL
	    e_el = c1 * c2 *r2;
#else
	    r = sqrt(r2);
	    e_el = c1 * c2 * r;
#endif
	    sum_e_el += e_el;
	}

	if (((a2 == 0)  && (b2 == 0)) || ((a1 == 0)  && (b1 == 0)))   {
	    e_lj6=0; e_lj12=0;
	}
	else {
	    r6 = r2 * r2 * r2 ; r12 = r6 * r6;
	    e_lj6 =   (a1 * a2 * r6);
	    e_lj12 =  (b1 * b2 * r12);
	    sum_e_lj += (e_lj12 - e_lj6);
	}
#undef DEBUGENER
#ifdef DEBUGENER
    r = sqrt(r2);
    fprintf(stderr,"%3d %s %3d %s %f %f %f\n",
    iat,pro[iat].atomName,jat,pro[jat].atomName, (-e_lj6+e_lj12)/4.182, e_el*331,1/r);
#endif
    }
    sum_e_el *= 331; sum_e_lj /= 4.182;
    ener = (REAL) (sum_e_el + sum_e_lj);

#ifdef SHOWPROGRESS
    fprintf(stderr," %7.2f %7.2f %7.2f kcal/mol\n",ener,sum_e_el,sum_e_lj);
#endif

    return (ener);
}

/******************************************************************
 *
 * AmidePairList() - Build list of atoms within NB_CUTOFF angstroms
 *                  from all specified atoms
 *
 ******************************************************************/

int AmidePairList(int *resorder_i)
{

int i,j,jat;
int numpairs=0;
int residue;

/* set the MOVE (=pro[].keydict) for the amide(s) in the group to 1 */
for (jat=0; jat<numAtms; jat++) MOVE(jat)=0;
i=0;
while ((residue=resorder_i[i])>-1) {
    residue++;
    for (jat=0; jat<numAtms; jat++) {
	if ((j=AMIDE(jat))!=0) 
	    if (j==residue || j==-residue) MOVE(jat)=j;
    }
    i++;
}

for (jat=0; jat<numAtms; jat++) {
    if (MOVE(jat)>0) AtomPairs(jat,&numpairs);
}

return (numpairs);

} /* end AmidePairList()   */

/******************************************************************
 *
 * AtomPairs() - Build list of atoms within NB_CUTOFF angstroms
 *                  from one specified atom
 *
 ******************************************************************/

void AtomPairs(int central_atom, int *numpairs)
{
int i,j,jat;
REAL adist,dist2;
int movei,movej;
static int allocate = 0;
static int new=TRUE;

if (new) {
    allocate=MAXPAIRS;
    NBPairs = (Pair *) malloc (MAXPAIRS * sizeof(Pair));
    if (!NBPairs) goto error;
    new = FALSE;
}
movei = MOVE(central_atom);

for(jat=0;jat<numAtms;jat++) {

  movej = MOVE(jat);
  if (movej) {
      if (movej<0) movej = -movej;
      if (movei==movej) continue; /* avoid self */
      if (jat > central_atom) continue; /* avoid duplicates */
  }
  else if (central_atom-jat>=1 && central_atom-jat<=2) continue; /* avoid CB-CD in gln and CA-CG in asn */

  /* OLD
  dist2=0.;
  for (j=0;j<3;j++)
      adist = pro[central_atom].xyz[j] - pro[jat].xyz[j];
      dist2 += adist*adist;
  }
  */

  dist2 = DistSq(pro[central_atom].xyz,pro[jat].xyz);
  
  if (dist2 > NB_CUTOFF_SQ)  continue; 
  
  NBPairs[*numpairs].iat = central_atom;
  NBPairs[*numpairs].jat = jat;
  (*numpairs)++;
  if (*numpairs >= allocate) {
     allocate *=2;
     NBPairs = (Pair *) realloc (NBPairs,allocate * sizeof (Pair));
     if (!NBPairs) goto error;
  }
}
 
 return ;

 error: fprintf(stderr,"* Failure to (re) allocate NBPairs\n"); exit (1);
 
}  /* end AtomPairs()   */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * int SetAmideFlag()
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int SetAmideFlag(PDB *pro_param)
{
int amNum, resnum;
int i,j,k,kk;
char *what, *wham;

static char *asn_atoms[]={
"CB", "CG", "OD1", "ND2", "HD1", "HD2", "" };
static char *gln_atoms[]={
"CG", "CD", "OE1", "NE2", "HE1", "HE2", "" };

char **amide_atoms;

  /* Get from Coordinate file all Asn and Gln Residue. Assign Consecutive Numbers 1 to amNum */
  resnum=-1;
  amNum = 0;
  for (i=0; i<numAtms;i++) AMIDE(i) = 0;
  for (i=0; i<numAtms;i++) {
    if (EQUAL(pro[i].resName,"ASN")) amide_atoms=asn_atoms;
    else if (EQUAL(pro[i].resName,"GLN")) amide_atoms=gln_atoms;
    else continue;

    if (pro[i].resSeq != resnum) {
	amNum ++;
	resnum = pro[i].resSeq;
    }

    what=pro[i].atomName;
    for (j=0; (TRUE) ;j++) {
	wham = amide_atoms[j];
	if (wham[0]==NULL_CHAR) break;
	if (EQUAL(wham,what)) break;
    }
    if (wham[0]==NULL_CHAR) continue; /* not part of the amide */
    if (j==0) AMIDE(i) = -amNum;
    else AMIDE(i) = amNum;
  }  /* end of loop over all protein atoms */

  /* save the atom indicess belonging to each amide */
  atoms_in_am = (int **) malloc (amNum * sizeof(int *));
  for (i=0;i<amNum;i++)  {
      atoms_in_am[i] = (int *) malloc (6 * sizeof(int));
      for (j=0;j<6;j++) atoms_in_am[i][j] = -1;
  }
  for (i=0;i<numAtms;i++) {
    if ((k=AMIDE(i))!=0) {
	if (k<0) k=-k; k--;
	for (kk=0;kk<6;kk++) 
	    if (atoms_in_am[k][kk]==-1) {
		atoms_in_am[k][kk]=i;
		break;
	    }
    }
  }

  return amNum;
} /* end SetAmideFlag() */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Search() search distancematrix for connectivities
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void Search(int residnum, int *prior)
{
    int i1, j1;

    /* Go through matrix by column */
    for (i1=0; i1<amNum; i1++) {
	if (distance[residnum][i1]==1) {
	    *prior=TRUE;
	    if (testcon(i1)==FALSE) { 
		connect[connum]=i1;
		connum++;
	    } 
	}
    }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * function check()
 *      check if residue is already on list
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int check(int number)
{
    int i1,j1;
    for (j1=0; j1<amNum; j1++) {
	for (i1=0; i1<amNum; i1++) {
	    if (resorder[j1][i1]==number) return 1;
	}
    }
    return 0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * testcon() Test if connectivity has already been dealt with
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int testcon(int rescon)
{
int j;
  for (j=0; j<=connum; j++) {
	 if (connect[j]==rescon) {
		return TRUE;
	 }
  }
  return FALSE;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * int extbit()
 *     Extract bits from bitstring:
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int extbit(int string,int bit)
{
int starter,bit2return;

  starter=1<<bit;
  bit2return=starter & string;
  bit2return=bit2return >> bit;
  return bit2return;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * int GroupAmides() 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int GroupAmides(void)
{
int i,j,res1,res2;
int amide1,amide2;
REAL atomdis;
int num, resnum;
int prior; /* Flag, TRUE when found connectivity */

  distance = (int **) malloc (amNum * sizeof(int *));
  for (i=0;i<amNum;i++)  {
      distance[i] = (int *) malloc (amNum * sizeof(int));
      if (!distance[i])  {
	  fprintf(stderr,"*** Unable to allocate memory for distance[%d]!\n",i);
	  exit(1); 
      }
  }
  for (i=0; i<amNum; i++) {
    for (j=0; j<amNum; j++) distance[i][j]=0;
  }

  /*  build Distance matrix */
  for (res1=0; res1<amNum; res1++) { 
    amide1 = atoms_in_am[res1][1];
    for (res2=res1+1; res2<amNum; res2++) {  
	amide2 = atoms_in_am[res2][1];
	atomdis = sqrt(DistSq(pro[amide1].xyz,pro[amide2].xyz));	
	if  (atomdis <= GROUP_CUTOFF)  distance[res1][res2]=distance[res2][res1]=1;
	else  distance[res1][res2]=0;
	/*		fprintf(stdout,"distance[%d][%d] = %f, (%d)\n",res1,res2,atomdis,distance[res1][res2]); */
	/*		fprintf(stdout,"amide1: %d, amide2: %d \n",amide1,amide2); */
    }

  }
#ifdef DEBUGX
    printf("Distancematrix: \n   ");
    for (i=0; i<amNum; i++) printf("%d ",i); printf("\n");
    for (i=0; i<amNum; i++) printf("..."); printf("\n"); 
    for (i=0; i<amNum; i++) {
	printf("%d .",i);
	for (j=0; j<amNum; j++) printf("%d ",distance[i][j]);
	printf("\n");
    } 

    printf("Sort Amides....\n");
#endif

  /* Sort Amides with respect to their proximity
	  in space.  Use array Resorder[res,resgroup_num] for sorting. */
  resorder = (int **) malloc ( (amNum + 1) * sizeof (int *));
  for (i=0;i< 2 * amNum * amNum;i++) {
    resorder[i]= (int *) malloc ( (amNum + 1) * sizeof (int ));
    if (!resorder[i])  {
	fprintf(stderr,"*** Unable to allocate memory for resorder[%d]!\n",i);
	exit(1);
    }
  }
  
  for (i=0;i< amNum + 1 ;i++) {
    for (j=0;j< amNum + 1;j++) resorder[i][j]=-1;
  }

  connum  = resgroup_num = num = 0; 
  for (resnum=0; resnum<amNum; resnum++ ) {
    Search(resnum, &prior);
    if (prior==TRUE) {
	prior=FALSE;
	if (check(resnum)==FALSE) {
#ifdef DEBUGX
	    fprintf(stdout,"Put amide %d in group %d as the %d s residue\n",resnum,resgroup_num,num);
#endif
	    resorder[resgroup_num][num]=resnum;
	    num++;
	}
	for (i=0; i<connum; i++) {
	    if (check(connect[i])==FALSE) { 			 
	    resorder[resgroup_num][num]=connect[i];
	    num++;
	}
	Search(connect[i], &prior);
    }
    resgroup_num++;
    num=0;
    connum=0;
    prior=FALSE;
#ifdef DEBUGX
    printf("Amide %d has been dealt with\n",resnum);
#endif
    }
#ifdef DEBUGX
    else printf("Amide %d has no connectivities \n",resnum);
#endif
  }

#ifdef DEBUGX
  printf("Deal now with nonconnected amides\n");
#endif
/* Put all remaining (nonconnected) residues on list */
  for (j=0; j<amNum; j++) {		
    if (check(j)==FALSE) {
	resorder[resgroup_num][num]=j;
	resgroup_num++;
#ifdef DEBUGX
	printf("Put amide %d as nonconnected on list\n",j);
#endif
    }
  }

#ifdef SHOWGROUPS
  fprintf(stdout,
      "* FLIP: Amides have been grouped according to proximity: \n");
  for (i=0; i<resgroup_num;i++) {
    if (resorder[i][0]==-1) continue;
    fprintf(stdout,"*  Amidegroup %d:\n",i);
    j=0;
    while ((res1=resorder[i][j])!=-1) {
	amide1 = atoms_in_am[res1][0];
	fprintf(stdout,"*      (%d) %d %s\n",res1,pro[amide1].resSeq,pro[amide1].resName);
	j++;
    }
  } 
#endif

return (resgroup_num);

} /* end GroupAmides() */
/* =========================================================================== */
