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
 *   FindHB program
 *
 *  Find potential hydrogen bonds.
 *  Input
 *      argv[1] - reformatted protein structure (output of reformatPDB)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define DISTSQ_MAX 12.25 /* 3.5 ^ 2 */
#define DISTSQ_MAXH 6.25 /* 2.5 ^ 2 */

#include "dowser.h"
#define HBTYPE(i) pro[i].key_dict
#define INDEX(i) pro[i].newResNo

extern void initVDW(char *filename, int num_atoms, PDB *atoms);
extern void readPDB10(char *filename, int *num_atoms, PDB **atoms_ptr);
extern REAL DistSq(REAL r1[3], REAL r2[3]);
int LocateMyPartners (int ihbtype, int ipro, int iindex);

PDB *pro;                       /* protein coordinates */
int numPro;
int ires;

void main(int argc, char *argv[])
{
int numPartners=0;
int index, hbtype, ipro;

if (argc != 2)  {
    fprintf(stderr,"correct usage: setOH input.pdb\n");
    exit(1);
}

/* read the protein coordinates */
readPDB10 (argv[1], &numPro, &pro);
initVDW(argv[1], numPro, pro);

/* loop over all the protein atoms to find all potential H-bonders */
for (ipro=0;ipro<numPro;ipro++)  {

    ires=pro[ipro].resSeq;
    index = -ires;
    hbtype = 0;
         if (EQUAL(pro[ipro].resName,"ser") && EQUAL(pro[ipro].atomName,"og"))  hbtype=3;
    else if (EQUAL(pro[ipro].resName,"thr")  && EQUAL(pro[ipro].atomName,"og1"))  hbtype=3;
    else if (EQUAL(pro[ipro].resName,"tyr")  && EQUAL(pro[ipro].atomName,"oh"))  hbtype=3;
    else if (EQUAL(pro[ipro].resName,"hoh")  && EQUAL(pro[ipro].atomName,"ow"))  hbtype=3;
    else if (EQUAL(pro[ipro].resName,"his")) {
	     if (EQUAL(pro[ipro].atomName,"n")) hbtype=1;
	else if (EQUAL(pro[ipro].atomName,"ce1")) { hbtype=3; index=ires; }
	else if (EQUAL(pro[ipro].atomName,"cd2")) { hbtype=3; index=ires; }
	else if (EQUAL(pro[ipro].atomName,"nd1")) { hbtype=3; index=ires; }
	else if (EQUAL(pro[ipro].atomName,"ne2")) { hbtype=3; index=ires; }
	else if (EQUAL(pro[ipro].atomName,"hd1")) { index=ires; }
    }
    else if (EQUAL(pro[ipro].atomType,"o")) hbtype=2;
    else if (EQUAL(pro[ipro].atomName,"n")) {
	if (EQUAL(pro[ipro].resName,"pro")) hbtype=0;
	else hbtype=1;
    }
    else if (EQUAL(pro[ipro].atomType,"n")) hbtype=1;
    else if (pro[ipro].LJ_c >= 1.) hbtype=1; /* positive ion */
    else if (pro[ipro].LJ_c <= -1.) hbtype=1; /* negative ion */
    
    HBTYPE(ipro) = hbtype;
    INDEX(ipro) = index;
}

/* for histidines only */
for (ipro=0;ipro<numPro;ipro++)  {
    if (EQUAL(pro[ipro].resName,"his") && HBTYPE(ipro)==3) {
	hbtype=HBTYPE(ipro);
	if (hbtype) numPartners +=  LocateMyPartners(hbtype,ipro,INDEX(ipro));
    }
}
fprintf (stdout,"Total number of partners found = %d\n",numPartners);

} /* end main() */


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *int LocateMyPartners (ihbtype,ipro,iindex)
 *  Find (potential) h-bonding partners for one atom
 *    The key to hbtype is as follows:
 *      hbtype = 0 : not a h-bonder
 *      hbtype = 1 : h-bond donor (no attached H)
 *      hbtype = 2 : h-bond acceptor (attached H has fixed position)
 *      hbtype = 3 : can be acceptor or donor, because attached H is movable.
 *    (Hydrogen atoms have hbtype=0)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int LocateMyPartners (int ihbtype, int ipro, int iindex)
{

int jpro;
REAL *ixyz = pro[ipro].xyz;
REAL *jxyz;
int count=0;
int jhbtype;
int foundone=FALSE;
REAL r2,r2h;
REAL *ihxyz, *jhxyz;
int use_hi, use_hj;
char flag;
int num_close;

/* find the position of the H atom */
if (ihbtype==1 && pro[ipro+2].back != ipro) { /* must be only 1 */
    ihxyz = pro[ipro+1].xyz; use_hi=TRUE;
}
else use_hi=FALSE;

num_close=0;
for (jpro=0;jpro<numPro;jpro++) {
    if (iindex==INDEX(jpro)) continue; /* avoid bonded atoms */
    if (!(jhbtype=HBTYPE(jpro))) { /* not h-bonder */
	if ((r2=DistSq(ixyz,pro[jpro].xyz)) < DISTSQ_MAX) num_close++;
	continue;
    }

    if (ihbtype==1 && jhbtype==1) continue; /* both donors */
    if (ihbtype==2 && jhbtype==2) continue; /* both acceptors */

    /* find the position of the H atom */
    if (jhbtype==1 && pro[jpro+2].back != jpro) { /* must be only 1 */
	jhxyz = pro[jpro+1].xyz; use_hj=TRUE;
    }
    else use_hj=FALSE;

    jxyz=pro[jpro].xyz;
    if ((r2=DistSq(ixyz,jxyz)) < DISTSQ_MAX) {
	flag = ' ';
	if (use_hi) {
	    if ((r2h=DistSq(ihxyz,jxyz)) > DISTSQ_MAXH) continue;
	    flag='H'; r2=r2h;
	}
	if (use_hj) {
	    if ((r2h=DistSq(jhxyz,ixyz)) > DISTSQ_MAXH) continue;
	    flag='H'; r2=r2h;
	}

	if (!foundone) {
	    fprintf(stdout,"For (%1d) %3d %3s %3s:",
		ihbtype,pro[ipro].resSeq,pro[ipro].resName,pro[ipro].atomName);
	    foundone=TRUE;
	}
	fprintf(stdout," (%1d) %5.1f %3d %3s %3s %c",
	    jhbtype,sqrt(r2),pro[jpro].resSeq,pro[jpro].resName,pro[jpro].atomName,
	    flag);
	count++;
    }
}

if (foundone) fprintf(stdout,"\n");
else {
    /*
    fprintf(stdout,"For %3d %3s %3s: %d near\n",
	pro[ipro].resSeq,pro[ipro].resName,pro[ipro].atomName,num_close);
    */
    
}

return (count);

} /* end LocateMyPartners() */
