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
 *   main: "dictmake"
 *      Program to create a dicionary for odd hetero residues
 *
 *   arguments: -pdbin FILENAME -resi RESNAME [-root root-atom-name]
 *   requires "dowser.h" and readPDB10()
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <math.h>
#include "dowser.h"

void readPDB10(char *filename, int *nAtom, PDB **atoms);

#define TRUE 1
#define MAXATOM 100
#define BONDMAXDIST 2.
#define MAXBOND 8
#define FILENAME "benzam.pdb"
#define NONSENSE -99.
#define PI180 57.295780
#define NOTSET -1

typedef struct {
    char name[5], resname[5], type[5]; int index, nbond, bond[MAXBOND]; REAL xyz[3];
    REAL length, angle, dihedral;
    int backchain, forwardchain;
} Atom;

void ReaRec();
char *NxtWord();
int  NxtInt();
REAL Distance(REAL *a, REAL *b);
REAL AngleValue(REAL *xi1, REAL *xi2, REAL *xi3);
REAL DihedralValue(REAL *xiat, REAL *xib, REAL *xibb, REAL *xibbb);
void CrossProduct(REAL *a, REAL *b, REAL *v);
REAL ScalarProduct3(REAL *a, REAL *b);
REAL Pair(REAL *a, REAL *b, REAL *u);
REAL DistSq(REAL *a, REAL *b);
void MakeChains(Atom *atoms, int nAtom);
void MakeBonds(Atom *atoms, int nAtom);
void CalcParam(Atom *atoms, int nAtom);
void ShowDict(Atom *atoms, int nAtom);
void DetermineType(Atom *atoms, int nAtom);

int main(int argc, char *argv[])
{
FILE *pdbfile;
char filename[256];
int eof;
int count=0, nAtom, nAtomPDB;
char *what;

Atom *atoms, *anatom, *anatomj;
Atom swapatom;
#define SWAP(x,y) swapatom=atoms[y]; for (iswap=y;iswap>x;iswap--) atoms[iswap]=atoms[iswap-1]; atoms[x]=swapatom;

PDB *PdbAtoms, *apdbatom;
int i,j,ib,nb,back,forw,jj,jjj,iswap;
int insertpoint;
REAL *xyz1, *xyz2;
REAL aaa;

int pdbin, resid, res_seq;
char rootatom[5]; int root;

    if (argc<2) {
	fprintf(stdout,"Normal usage: dictmake -pdbin pdbfile -residue residue-name [-root name-of-root-atom]\n");
	exit (1);
    }
    /* look for arguments that affect the way dowser is to be used */
    *rootatom = NULLCHAR;
    for (i=1; i<argc; i++) {
	if (EQUAL(argv[i],"-root")) strncpy(rootatom,argv[i+1],4);
	if (EQUAL(argv[i],"-pdbin")) pdbin=i+1;
	if (EQU(argv[i],"-resi")) resid=i+1;
	/* if (EQUAL(argv[i],"-pdbout")) pdbout=i+1; */
    }

    /* process the input PDB file, which has been grepped to get only one new molecule  */
    /* the readPDB routine allocates the needed storage for atoms[] */
    readPDB10(argv[pdbin], &nAtomPDB, &PdbAtoms);

    /* Count the occurrence of the selected residue in the pdb file */
    res_seq=0;
    nAtom=0;
    apdbatom=PdbAtoms;
    for (i=0;i<nAtomPDB;i++) {
	if (EQUAL(apdbatom->resName,argv[resid])) {
	    /* proces only one residue of the new type */
	    if (res_seq==0) res_seq = apdbatom->resSeq;
	    else if (res_seq != apdbatom->resSeq) break; 
	    nAtom++;
	}
	apdbatom++;
    }
    fprintf (stdout, "REMARK %d atoms of type %s\n",nAtom,argv[resid]);
    if (nAtom==0) exit (1);

    atoms = (Atom *) malloc(nAtom*sizeof(Atom));

    /* extract the information from the pdb records */
    apdbatom=PdbAtoms;
    anatom=atoms;
    for (i=0;i<nAtomPDB;i++) {
	if (EQUAL(apdbatom->resName,argv[resid])) {
	    strcpy(anatom->name,apdbatom->atomName);
	    strcpy(anatom->resname,apdbatom->resName);
	    anatom->xyz[0]=apdbatom->xyz[0];
	    anatom->xyz[1]=apdbatom->xyz[1];
	    anatom->xyz[2]=apdbatom->xyz[2];
	    anatom->nbond=0;
	    anatom->length=anatom->angle=anatom->dihedral=0.;
	    anatom->backchain=anatom->forwardchain=NOTSET;
	    for (j=0;j<MAXBOND;j++) anatom->bond[j]=-99;
	    anatom++;
	}
	apdbatom++;
    }

/* ======================================================================= */
    MakeBonds(atoms,nAtom);
/* ======================================================================= */
    MakeChains(atoms,nAtom);
/* ======================================================================= */
#undef TEST
#ifdef TEST
    CalcParam(atoms,nAtom);
    ShowDict(atoms,nAtom);
#endif
/* ======================================================================= */
    /* Reorder the residue */
    if (*rootatom) {
	/* place the root atom in position 0 */
	root = NOTSET;
	for (i=0;i<nAtom;i++) {
	    if (EQUAL(atoms[i].name,rootatom)) { root = i; }
	    atoms[i].index = i;
	}
	if (root == NOTSET) { 
	    fprintf(stdout,"The requested root atom %s is not present\n",rootatom);
	    exit (1);
	}
	if (root !=0 ) {
	    SWAP(0,root);
#undef SHOWSWAP
#ifdef SHOWSWAP
	    fprintf (stdout,"SWAP(%d,%d)\n",0,root);
#endif
	}
    }
    else root=0;

    /* order the remaining atoms */
    for (i=0;i<nAtom;i++) {
	insertpoint = i+1;
	anatom = atoms + i;
	back = anatom->index; 
	for (j=i+1;j<nAtom;j++) {
	    if (j<insertpoint) continue;
	    anatomj=atoms+j;
	    nb = anatomj->nbond;
	    for (ib=0;ib<nb;ib++) {
		if (back == anatomj->bond[ib]) {
		    if (insertpoint!=j) {
			SWAP (insertpoint,j);
#ifdef SHOWSWAP
			fprintf (stdout,"SWAP(%d,%d)\n",insertpoint,j);
#endif
		    }
		    insertpoint++;
		    break;
		}
	    }
	}

    }

    MakeBonds(atoms,nAtom);
    MakeChains(atoms,nAtom);
/* ====================================================================== */

    CalcParam(atoms,nAtom);
    DetermineType(atoms,nAtom);
    ShowDict(atoms,nAtom);

} /* end main() */

/*
REMARK ATOM resn atom back forward bond angle dihedral charge type
RESIDUE ALA TERM NH3 COO
ATOM ALA  N    C    CA    1.320  114.0  180.0  -0.280 N
ATOM ALA  H    N    NOT   1.000  123.0    0.0   0.280 H
*/

/* * * * * * * * * * * * * * * * * * *
 *  REAL Distance
 * * * * * * * * * * * * * * * * * * */
REAL Distance (REAL *a, REAL *b)
{
    int i;
    REAL aaa=0.;
    REAL bbb;

    for (i=0;i<3;i++) {
	bbb = a[i] - b[i];
	aaa += bbb*bbb;
    }
    return sqrt(aaa);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  REAL DihedralValue(x,iat,ib,ibb,ibbb)
 *      x= coordinate set
 *	iat,ib,ibb,ibbb = four atom indices, (start at 1)
 *	returns dihedral angle in radians (each "b" is a backchain)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
REAL DihedralValue(REAL *xiat,REAL *xib,REAL *xibb,REAL *xibbb)
{
REAL v[3],w[3],u[3],e1[3],e2[3],ee[3];
REAL wl,vl,ul,cr,sr;

    wl=Pair(xib,xiat,w);
    vl=Pair(xibb,xib,v) ;
    ul=Pair(xibbb,xibb,u) ;
    CrossProduct(v,w,e2) ;
    CrossProduct(u,v,e1);
    cr=ScalarProduct3(e1,e2) ; CrossProduct(e1,e2,ee) ;
    sr=ScalarProduct3(ee,v)/vl;
    return (atan2(sr,cr));
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  REAL AngleValue(x,i1,i2,i3)
 *      x= coordinate set
 *	i1,i2,i3 = three atom indices, (start at 1)
 *	returns angle i1-i2-i3, in radians
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
REAL AngleValue(REAL *xi1,REAL *xi2,REAL *xi3)
{
REAL d1,d2,d3;
d1=Distance(xi1,xi2);
d2=Distance(xi2,xi3);
d3=Distance(xi1,xi3);
return(acos(-(d3*d3-d1*d1-d2*d2)/(2.*d1*d2)) );
}

/* * * * * * * * * * * * * * * * * * * * * * * * *
 *  CrossProduct(a,b,v)
 *     Vector cross product (v = a x b)
 * * * * * * * * * * * * * * * * * * * * * * * * */
void CrossProduct(REAL *a,REAL *b,REAL *v)
 {
     int i,j,k;
     for (i=0;i<3;i++) {
         j=i+1 ; if (j==3) j=0; k=3-i-j;
         *(v+i)= *(a+j)* *(b+k) - *(a+k)* *(b+j);
     }
 }


/* * * * * * * * * * * * * * * * * * * * * * * * *
 * REAL ScalarProduct3(a,b)
 *     Returns vector scalar product (a . b)
 * * * * * * * * * * * * * * * * * * * * * * * * */
 REAL ScalarProduct3(REAL *a,REAL *b)
 {
     return(*a* *b + *(a+1)* *(b+1) + *(a+2)* *(b+2));
 }

/* * * * * * * * * * * * * * * * * * * * * * * * * *
 *   dist_sq: returns distance between 2 points
 * * * * * * * * * * * * * * * * * * * * * * * * * */
REAL DistSq(REAL *a,REAL *b)
{
    float r2;
    int k;

    r2=0.;
    for (k=0;k<3;k++)
        r2+=(*(a+k)-*(b+k))*(*(a+k)-*(b+k));
    return(r2);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * REAL Pair(a,b,u)
 *      a=coordinates of point 1
 *      b=coordinates of point 2
 *      u=vector b-a
 *      function returns the distance between a and b
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
REAL Pair(REAL *a,REAL *b,REAL *u)
{
int i;
/*
  a,b = 2 positions (input), u=vector a->b (output), pair = |u|
*/
for (i=0;i<3;i++) u[i]=b[i]-a[i];
return (sqrt(ScalarProduct3(u,u)));
}

/* * * * * * * * * * * * * * * * * * *
 *  MakeChains()
 * * * * * * * * * * * * * * * * * * */

void MakeChains(Atom *atoms, int nAtom)
{
int i,j;
Atom *anatom, *anatomj;
int forw, back, ib, nb;

    /* Determine the Backchains */
    anatom=atoms;
    for (i=0;i<nAtom;i++) {
	forw = back = NOTSET;
	/* loop over all prior atoms */
	for (j=0;j<i;j++) {
	    anatomj = atoms+j;
	    nb = anatomj->nbond;
	    for (ib=0;ib<nb;ib++) {
		if (anatomj->bond[ib]==i) { back=j; goto FoundIt; }
	    }
	}
	FoundIt: anatom->backchain=back;
	anatom++;
    }

    /* Determine the Forwardchains */
    anatom=atoms;
    for (i=0;i<nAtom;i++) {
	/* loop over all subsequent atoms */
	/* highest forward chain */
	forw=NOTSET;
	for (j=i+1;j<nAtom;j++) {
	    anatomj = atoms+j;
	    if (anatomj->backchain==i) forw=j; 
	}
	anatom->forwardchain=forw;
	anatom++;
    }
    return;
}

/* * * * * * * * * * * * * * * * * * *
 *  MakeBonds()
 * * * * * * * * * * * * * * * * * * */
void MakeBonds(Atom *atoms, int nAtom)
{
int i,j;
Atom *anatom, *anatomj;
int forw, back, ib, nb;
REAL *xyz1, *xyz2;
REAL aaa;

    anatom=atoms;
    for (i=0;i<nAtom;i++) {
	anatom->nbond=0;
	anatom->length=anatom->angle=anatom->dihedral=0.;
	anatom->backchain=anatom->forwardchain=NOTSET;
	for (j=0;j<MAXBOND;j++) anatom->bond[j]=-99;
	anatom++;
    }
    /* Find the Bonds */
    for (i=0;i<nAtom;i++)  {
	xyz1=&(atoms[i].xyz[0]);
	for (j=i+1;j<nAtom;j++) {
	    xyz2=&(atoms[j].xyz[0]);
	    aaa = Distance (xyz1,xyz2);
	    if (aaa < BONDMAXDIST) {
#ifdef TEST
		fprintf (stdout,"Bond %s-%s %8.2f\n",atoms[i].name,atoms[j].name,aaa);
#endif
		if (atoms[i].nbond == MAXBOND) {
		   fprintf(stdout,"Too many bonds to one atom\n"); exit(1);
		}
		if (atoms[j].nbond == MAXBOND) {
		   fprintf(stdout,"Too many bonds to one atom\n"); exit(1);
		}
		atoms[i].bond[atoms[i].nbond] = j;
		atoms[j].bond[atoms[j].nbond] = i;
		++(atoms[i].nbond);
		++(atoms[j].nbond);
	    }
	}

    }
#undef SHOW1
#ifdef SHOW1
    anatom=atoms;
    fprintf(stdout,"RESIDUE %-4s\n",anatom->resname);
    for (i=0;i<nAtom;i++) {
       fprintf(stdout,"% 2d %-4s %8.3f %8.3f %8.3f",
	   i,atoms[i].name,atoms[i].xyz[0],atoms[i].xyz[1],atoms[i].xyz[2]);
	   for (j=0;j<atoms[i].nbond;j++) fprintf(stdout," ->%-4s",atoms[atoms[i].bond[j]].name);
	   fprintf(stdout,"\n");
    }
#endif
}

/* * * * * * * * * * * * * * * * * * *
 *  ShowDict()
 * * * * * * * * * * * * * * * * * * */
void ShowDict(Atom *atoms, int nAtom)
{
int i;
Atom *anatom;

    fprintf (stdout,
	"REMARK resn atom back forward bond angle dihedral charge type\n");
    anatom=atoms;
    fprintf(stdout,"RESIDUE %-4s\n",anatom->resname);
    for (i=0;i<nAtom;i++) {

       fprintf(stdout,"ATOM %-4s %-4s %-4s %-4s %8.3f %8.1f %8.1f",
	   anatom->resname,
	   anatom->name,
	   (anatom->backchain>NOTSET)? atoms[anatom->backchain].name: "NOT",
	   (anatom->forwardchain>NOTSET)? atoms[anatom->forwardchain].name: "NOT",
	   anatom->length,
	   (anatom->angle!=NONSENSE)?PI180*anatom->angle:anatom->angle,
	   (anatom->dihedral!=NONSENSE)?PI180*anatom->dihedral:anatom->dihedral);

	   fprintf(stdout,"  ?.??? %-s",anatom->type);
	   fprintf(stdout,"\n");
	anatom++;
    }
}

/* * * * * * * * * * * * * * * * * * *
 *  CalcParam()
 * * * * * * * * * * * * * * * * * * */
void CalcParam(Atom *atoms, int nAtom)
{
int i,j,jj,jjj;
Atom *anatom;
    /* find the geometric parameters */
    anatom = atoms;
    for (i=0;i<nAtom;i++) {
	anatom->length = anatom->angle = anatom->dihedral = NONSENSE;
	j = anatom->backchain;
	if (j>NOTSET) {
	     anatom->length = Distance(anatom->xyz,atoms[j].xyz);
	}
	else { anatom++; continue; }
	jj = atoms[j].backchain;
	if (jj>NOTSET) {
	     anatom->angle = AngleValue(anatom->xyz,atoms[j].xyz,atoms[jj].xyz);
	}
	else { anatom++; continue; }
	jjj = atoms[jj].backchain;
	if (jjj>NOTSET) {
	     anatom->dihedral =
		 DihedralValue(anatom->xyz,atoms[j].xyz,atoms[jj].xyz,atoms[jjj].xyz);
	}
	anatom++;
    }
}

/* * * * * * * * * * * * * * * * * * *
 *  DetermineType()
 *    determine type of carbon atoms on the basis of geometric criteria
 * * * * * * * * * * * * * * * * * * */
void DetermineType(Atom *atoms, int nAtom)
{
int i,j,jj,jjj, j1, j2;
int n_change;
Atom *anatom, *anatomj, *anatomk, *anatom1, *anatom2;
char *atomtype;
REAL aaa;
int back,forw,link;

    /* intialize type to one character plus '?' */ 
    anatom=atoms;
    for (i=0;i<nAtom;i++) {
	atomtype=anatom->type;
	atomtype[0]=anatom->name[0];
	if (*atomtype == 'C' || *atomtype == 'O') {
	    atomtype[1] = '?';
	    atomtype[2] = NULLCHAR;
	}
	else atomtype[1]=NULLCHAR;
	anatom++;
    }

    /* look for dihedrals that are zero, which indicates aromaticity */
    anatom=atoms;
    for (i=0;i<nAtom;i++) {
	if (anatom->dihedral>-10./PI180 && anatom->dihedral<10./PI180) {
	    anatomj = atoms + anatom->backchain;
	    anatomk = atoms + anatomj->backchain;

	    if (*(anatomj->type)=='C') {
		if (anatomj->nbond==3) {
		    strcpy(anatomj->type,"CR");
		}
		else if (anatomj->nbond==2) {
		    strcpy(anatomj->type,"CHR");
		}
	    }
	    if (*(anatomk->type)=='C') {
		if (anatomk->nbond==3) {
		    strcpy(anatomk->type,"CR");
		}
		else if (anatomk->nbond==2) {
		    strcpy(anatomk->type,"CHR");
		}
	    }
	}
	anatom++;
    }

    /* look at value of bondangles which may indicate aromaticity */
    anatom=atoms;
    for (i=0;i<nAtom;i++) {
	if (anatom->angle>115./PI180 && anatom->angle<125./PI180) {
	    anatomj = atoms + anatom->backchain;

	    if (EQUAL(anatomj->type,"C?")) {
		if (anatomj->nbond==3) {
		    strcpy(anatomj->type,"CR"); continue;
		}
		else if (anatomj->nbond==2) {
		    strcpy(anatomj->type,"CHR"); continue;
		}
	    }
	}
	anatom++;
    }

    /* Look for broken chain, indicates ring, may be aromatic */
    anatom=atoms;
    for (i=0;i<nAtom;i++) {
	if (EQUAL(anatom->type,"C?") && anatom->nbond > 1) {
	    if ((back=anatom->backchain)==NOTSET || (forw=anatom->forwardchain)==NOTSET) {
		if (back==NOTSET) jj=forw; else jj=back;
		j=0;
		while (j<anatom->nbond && (link=anatom->bond[j])!=jj) {
		    anatomj=atoms+link;
		    for (j1=0;j1<anatom->nbond;j1++) {
			anatom1=atoms+anatom->bond[j1];
			for (j2=0;j2<anatomj->nbond;j2++) {
			    anatom2=atoms+anatomj->bond[j2];
			    aaa=DihedralValue(anatom->xyz,anatom1->xyz,anatom2->xyz,anatomj->xyz);
			    if (aaa > -10./PI180 && aaa < 10./PI180) {
				if (anatom->nbond==3) {
				    strcpy(anatom->type,"CR"); goto NextAtom;
				}
				else if (anatom->nbond==2) {
				    strcpy(anatom->type,"CHR"); goto NextAtom;
				}
			    }
			}
		    }
		    j++;
		}
	    }
	}
	NextAtom: anatom++;
    }

    /* look for planar group: one carbon atom plus three substituents */
    for (i=0;i<nAtom;i++) {
	anatom = atoms +i;
	if (EQUAL(anatom->type,"C?") && anatom->nbond == 3) {
	    anatomj=atoms + anatom->bond[0];
	    anatom1=atoms + anatom->bond[1];
	    anatom2=atoms + anatom->bond[2];
	    aaa=DihedralValue(anatomj->xyz,anatom->xyz,anatom1->xyz,anatom2->xyz);
	    if (aaa > PI-10./PI180 || aaa < -PI+10./PI180) 
		strcpy(anatom->type,"CR");
	}
    }

    /* process the remaining C atoms without type as aliphatic */
    anatom=atoms;
    for (i=0;i<nAtom;i++) {
	if (EQUAL(anatom->type,"C?")) {
	    if (anatom->nbond==4) {
		strcpy(anatom->type,"CA");
	    }
	    else if (anatom->nbond==3) {
		strcpy(anatom->type,"CH1");
	    }
	    else if (anatom->nbond==2) {
		strcpy(anatom->type,"CH2");
	    }
	}
	anatom++;
    }

}
/*
REMARK atomtype LJ-a LJ-b    MS-radius
TYPE   N      49.36   1300.0  2.40   # any Nitorgen
TYPE   H       0.00      0.0  0.0    # polar hydrogen
TYPE   CR     48.37   1837.0  2.89   # sp2 carbon, w/o H
TYPE   OA     47.56   1230.0  2.20   # oxygen in carboxylate
TYPE   O      47.56    861.0  2.20   # other oxygen
TYPE   S      99.20   3617.0  2.70   # any Sulfur
TYPE   OW     51.16   1623.0  2.20   # SPC water oxygen
REMARK: ions
TYPE   CAL    47.56   861.00  1.80   # calcium
TYPE   NA     47.56   861.00  1.80   # sodium
TYPE   ZN     20.57   942.38  1.60   # zinc
REMARK: united-atom parameters for C and H  
TYPE   CH2    95.38   5944.0  3.00   # aliphatic C with 2 H
TYPE   CH1   111.80   8470.4  3.00   # aliphatic C with 1 H
TYPE   CH3    94.20   5114.0  3.00   # aliphatic C with 3 H
TYPE   CHR    74.25   3888.0  2.89   # aromatic C with 1 H
REMARK: all-atom parameters for C and H  
TYPE   HA      9.20   123.0   1.90   # aliphatic H
TYPE   CA     40.0   1959.0   2.90   # aliphatic C
*/
