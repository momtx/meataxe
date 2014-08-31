/* ============================= C MeatAxe ==================================
   File:        $Id: mksub.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Calculate the submodule lattice.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <string.h>
#include <stdlib.h>




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

#define O_MOUNTAINS		0x01
#define O_SUBMODULES		0x02
#define O_DOTTEDLINES		0x04
#define O_EXTFILES		0x08
#define O_RADICAL		0x10
#define O_SOCLE			0x20
#define O_INCIDENCES		0x40
#define O_ALL (O_MOUNTAINS|O_SUBMODULES|O_DOTTEDLINES|O_EXTFILES|\
	       O_RADICAL|O_SOCLE|O_INCIDENCES)


int opt_b = 0;			/* -b option (blocks) */
int opt_o = O_ALL;
int opt_G = 0;
int done[LAT_MAXCF];		/* */
int blnum;			/* Number of current block */
int blsize;			/* Block size */
int block[LAT_MAXCF];		/* Block members */

int firstm[LAT_MAXCF+1];		/* First mountain */
int firstdl[LAT_MAXCF+1];		/* First dotted line */


/* Data read from input files
   -------------------------- */
int xnmount = 0;		/* Number of mountains */
int xndotl = 0;			/* Number of dotted lines */
BitString_t *xsubof[MAXCYCL];	/* Incidence matrix */
BitString_t *xdotl[MAXDOTL];	/* Dotted lines */
long xmdim[MAXCYCL];		/* Mountain dimensions */
static Lat_Info LI;		/* Data from .cfinfo file */


/* Data for current block
   ---------------------- */
int bnmount, bndotl;		/* As above, but for current block */
BitString_t *bsubof[MAXCYCL];
BitString_t *bsupof[MAXCYCL];	/* Transposed incidence matrix */
BitString_t *bdotl[MAXDOTL];
BitString_t *bdlspan[MAXDOTL];
long bmdim[MAXCYCL];

/* Data used during computation
   ---------------------------- */
BitString_t *sub[MAXNSUB];		/* Submodules */
int nsub = 0;			/* Number of submodules */
int lastgen = 0;		/* First submodule of last generation */
int generation;			/* Generation count */
int nadd;			/* Number of calls to addtolist() */
int oldnsub;			/* */
BitString_t *y;			/* Temporary bit string */
long *subdim;			/* Submodule dimensions */
char *israd;			/* Radical series */
char *issoc;			/* Socle series */
char *ismount;			/* Mountains */
int **max;			/* List of maximal submodules */



static MtxApplicationInfo_t AppInfo = { 
"mksub", "Find Submodules",
"\n"
"SYNTAX\n"
"    mksub [<Options>] <Name>\n"
"\n"
"ARGUMENTS\n"
"    <Name> .................. Name of the representation\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... GAP output (implies -Q)\n"
"    -b ...................... Find blocks\n"
"    -o <Fmt> ................ Output elements in <Fmt>\n"
"    -n <Fmt> ................ Exclude elements in <Fmt>\n"
"        <Fmt> is any combination of m (mountains), d (dotted lines),\n"
"        i (incidence matrix), e (.lat and .gra files), s (submodule list),\n"
"        r (radical series), and o (socle series).\n"
"\n"
"FILES\n"
"    <Name>.cfinfo ........... IO Constituent info file\n"
"    <Name>.inc .............. I  Incidence matrix generated by MKINC\n"
"    <Name>.dot .............. I  Dotted-lines generated by MKDOTL\n"
"    <Name>.mnt .............. I  Mountain dimensions\n"
"    <Name>.out .............. O  Submodule lattice\n"
"    <Name>.lat .............. O  Incidence matrix of the submodules (GAP)\n"
"    <Name>.gra .............. O  Submodule lattice for MKGRAPH\n"
"\n"
"    If -b is used, output files are produced for each block, and a\n"
"    block number is appended to the file names (e.g., `psl27.out.1').\n"
};


static MtxApplication_t *App = NULL;




/* -----------------------------------------------------------------
   init() - Read .cfinfo and .inc file
   ----------------------------------------------------------------- */

static void init(const char *basename)

{	
    int i;
    long l[2];
    char fn[40];
    FILE *f;

    if (Lat_ReadInfo(&LI,basename) != 0)
    {
	MTX_ERROR1("Error reading %s.cfinfo",basename);
	return;
    }
    
    /* Read incidence matrix
       --------------------- */
    f = SysFopen(strcat(strcpy(fn,LI.BaseName),".inc"),FM_READ);
    if (f == NULL) 
    {
	MTX_ERROR1("Cannot open %s",fn);
	return;
    }
    SysReadLong(f,l,1);
    xnmount = (int) l[0];
    MESSAGE(1,("Reading%s: %d mountain%s\n",fn,xnmount,xnmount == 1 ? "" : "s"));
    if (xnmount > MAXCYCL) 
    {
	MTX_ERROR2("Too many mountains (%d, max=%d)",xnmount,MAXCYCL);
	return;
    }
    for (i = 0; i < xnmount; ++i)
    {
	if ((xsubof[i] = BsRead(f)) == NULL)
	{
	    MTX_ERROR("Error reading incidence matrix");
	    return;
	}
	if (xsubof[i]->Size != xnmount)
	{
	    MTX_ERROR("Invalid bit string");
	    return;
	}
    }
    fclose(f);

    /* Read dotted lines
       ----------------- */
    f = SysFopen(strcat(strcpy(fn,LI.BaseName),".dot"),FM_READ);
    if (f == NULL) 
    {
	MTX_ERROR1("Cannot open %s",fn);
	return;
    }
    MESSAGE(1,("Reading %s: ",fn));
    SysReadLong(f,l,1);
    xndotl = (int) l[0];
    MESSAGE(1,("%d dotted line%s\n",xndotl,xndotl == 1 ? "" : "s"));
    if (xndotl > MAXDOTL) 
    {
	MTX_ERROR2("Too many dotted-lines (%d, max=%d)",xndotl,MAXDOTL);
	return;
    }
    for (i = 0; i < xndotl; ++i)
    {
	if ((xdotl[i] = BsRead(f)) == NULL)
	{
	    MTX_ERROR("Error reading dotted lines");
	    return;
	}
    }
    fclose(f);
    y = BsAlloc(xnmount);

    /* Read dimensions
       --------------- */
    f = SysFopen(strcat(strcpy(fn,LI.BaseName),".mnt"),FM_READ|FM_TEXT);
    if (f == NULL) 
    {
	MTX_ERROR1("Cannot open %s",fn);
	return;
    }
    MESSAGE(1,("Reading %s\n",fn));
    for (i = 0; i < xnmount; ++i)
    {
	long mno, mdim;
	if (fscanf(f,"%ld%ld",&mno,&mdim) != 2 || mno != i || mdim < 1)
	{
	    MTX_ERROR("Error in .mnt file");
	    return;
	}
	xmdim[i] = mdim;
	while (getc(f) != '\n')	/* Skip class */
	{
	    if (ferror(f) || feof(f))
	    {
		MTX_ERROR("Error in .mnt file");
		return;
	    }
	}
    }
    fclose(f);
}


/* -----------------------------------------------------------------
   init2()
   ----------------------------------------------------------------- */

static void init2()

{
    int i;

    /* Set firstm and firstdl
       ---------------------- */
    firstm[0] = 0;
    firstdl[0] = 0;
    for (i = 0; i < LI.NCf; ++i)
    {
	firstm[i+1] = firstm[i] + LI.Cf[i].nmount;
	firstdl[i+1] = firstdl[i] + LI.Cf[i].ndotl;
    }

    /* Initialize done[]
       ----------------- */
    for (i = 0; i < LI.NCf; ++i) done[i] = 0;
}


/* -----------------------------------------------------------------
   isotype() - Finde den Isomorphietyp eines gegebenen Mountains.
   ----------------------------------------------------------------- */

static int isotype(int mnt)

{
    int m;

    for (m = 0; (mnt -= LI.Cf[block[m]].nmount) >= 0; ++m);
    return (block[m]);
}


/* -----------------------------------------------------------------
   findrsm() - Find the radical and socle series, and the
	mountains.
   ----------------------------------------------------------------- */

void findrsm()

{
    char *flag = NALLOC(char,nsub);
    BitString_t *bs = BsAlloc(bnmount);
    int i, k;


    /* Berechne maximale Teilmoduln und Dimensionen
       -------------------------------------------- */
    ismount = NALLOC(char,nsub);
    max = NALLOC(int *,nsub);
    subdim = NALLOC(long,nsub);

    for (i = 0; i < nsub; ++i)
    {
	int maxcount = 0, *lp;
        memset(flag,0,(size_t) nsub);

	/* Bestimme alle maximalen Teilmoduln
	   ---------------------------------- */
	for (k = i-1; k >= 0; --k)
	{
	    if (flag[k] != 0) continue;
	    if (BsIsSub(sub[k],sub[i]))
	    {
		int l;
		flag[k] = 1;
		++maxcount;
		for (l = k-1; l >= 0; --l)
		{
		    if (BsIsSub(sub[l],sub[k])) 
			flag[l] = 2;
		}
	    }
	}

	/* Baue eine Liste mit den Nummern der maximalen Teilmoduln
	   und den Isomorphietypen der einfachen Faktoren auf.
	   -------------------------------------------------------- */
	lp = max[i] = NALLOC(int,2*maxcount+1);
	for (k = 0; k < i; ++k)
	{
	    if (flag[k] == 1)
	    {
		int l;
		*lp++ = k;
		for (l = 0; !BsTest(sub[i],l)||BsTest(sub[k],l); ++l);
		*lp++ = isotype(l);
	    }
	}
	*lp = -1;
	ismount[i] = (char) (maxcount == 1);

	/* Berechne die Dimension. Da wir aufsteigend vorgehen,
	   ist die Dimension der Maximalen hier schon bekannt.
	   ---------------------------------------------------- */
	if (maxcount == 0)
	    subdim[i] = 0;
	else
	{
	    subdim[i] = subdim[max[i][0]] + LI.Cf[max[i][1]].dim;
	}
    }

    /* Berechne die Radikalreihe
       -------------------------- */
    israd = NALLOC(char,nsub);
    memset(israd,0,(size_t)nsub);
    for (i = nsub-1; i > 0; )
    {	
	int *lp;
	BsCopy(bs,sub[i]);
	for (lp = max[i]; *lp >= 0; lp += 2)
	    BsAnd(bs,sub[*lp]);
	for (i = nsub-1; !BsIsSub(sub[i],bs); --i);
	israd[i] = 1;
    }

    /* Berechne die Sockelreihe
       ------------------------ */
    issoc = NALLOC(char,nsub);
    memset(issoc,0,(size_t)nsub);
    for (i = 0; i < nsub-1; )
    {
	/* Find the simple submodules
	   -------------------------- */
	memset(flag,0,(size_t) nsub);
	for (k = i+1; k < nsub; ++k)
	{
	    if (flag[k] != 0) continue;
	    if (BsIsSub(sub[i],sub[k]))
	    {
		int l;
		flag[k] = 1;
		for (l = k + 1; l < nsub; ++l)
		    if (BsIsSub(sub[k],sub[l])) flag[l] = 2;
	    }
        }

        /* Calculate the sum of all simple submodules (=Socle)
	   --------------------------------------------------- */
	BsCopy(bs,sub[i]);
	for (k = i; k < nsub; ++k)
	{
	    if (flag[k] == 1) 
		BsOr(bs,sub[k]);
	}

	for (i = 0; !BsIsSub(bs,sub[i]); ++i);
	issoc[i] = 1;
    }
}


/* -----------------------------------------------------------------
   extend() - Add one mountain to a given BitString_t
   ----------------------------------------------------------------- */

static char dlflag[MAXDOTL] = {0};

static void extend(BitString_t *x, int i, int nextend)

{
    int k;
    int changed;

    BsOr(x,bsupof[i]);		/* Add mountain and its subspaces */
    if (nextend) BsClear(x,i);	/* Add the radical only */

    /* Make closure
       ------------ */
    memset(dlflag,0,sizeof(dlflag));
    for (changed = 1; changed; )
    {
	changed = 0;
        for (k = 0; k < bndotl; ++k)
        {
	    if (!dlflag[k] && BsIntersectionCount(x,bdotl[k]) >= 2)
	    {
	        BsOr(x,bdlspan[k]);
	        dlflag[k] = 1;
	        changed = 1;
	    }
	}
    }
}


/* -----------------------------------------------------------------
   writeresult() - Write output files
   ----------------------------------------------------------------- */

FILE *openout(char *name)

{
    FILE *f;
    char fn[200];

    sprintf(fn,opt_b ? "%s%s.%d" : "%s%s",LI.BaseName,name,blnum);
    MESSAGE(1,("Writing %s\n",fn));
    f = SysFopen(fn,FM_CREATE|FM_TEXT);
    if (f == NULL)
    {
	MTX_ERROR1("Cannot open %s",fn);
	return NULL;
    }
    return f;
}

#define NDIG(x) (x>99999?6:x>9999?5:x>999?4:x>99?3:x>9?2:1)

static int printbs(FILE *f, BitString_t *b)

{
    int k, k1, k2, len, flag;

    if (bnmount < 100)
    {
	for (k = 0; k < bnmount;  ++k)
	    putc(BsTest(b,k) ? '+' : '.',f);
	len = bnmount;
    }
    else
    {
	len = k = 0;
	flag = 0;
	while (1)
	{
	    while (k < bnmount && !BsTest(b,k)) ++k;
	    if (k >= bnmount) 
		break;
	    k1 = k;
	    while (k < bnmount && BsTest(b,k)) ++k;
	    k2 = k-1;
	    if (flag)
	    {
		putc(',',f);
		++len;
	    }
	    else
		flag = 1;
	    if (k2 > k1)
	    {
		fprintf(f,"%d-%d",k1,k2);
		len += NDIG(k1)+NDIG(k2)+1;
	    }
	    else
	    {
		fprintf(f,"%d",k1);
		len += NDIG(k1);
	    }
	}
    }
    return len;
}



void writeresult()

{
    FILE *f, *g;
    int i, k;
    char tmp[100];
    BitString_t *b = BsAlloc(bnmount);

    MESSAGE(0,("Finished, %d submodules found\n",nsub));
    f = openout(".out");

    /* Write irreducibles
       ------------------ */
    fprintf(f,"Irreducibles:\n");
    fprintf(f,"    Type   Mult   SF   Mountains           Dotted lines\n");
    for (i = 0; i < blsize; ++i)
    {
	sprintf(tmp,"%s",Lat_CfName(&LI,block[i]));
	fprintf(f,"    %-7s%-7ld%-5ld",tmp,LI.Cf[block[i]].mult,
	    LI.Cf[block[i]].spl);
	sprintf(tmp,"%ld (%ld-%ld)",LI.Cf[block[i]].nmount,
	    (long) firstm[block[i]],
	    LI.Cf[block[i]].nmount+firstm[block[i]]-1);
	fprintf(f,"%-20s",tmp);
	if (LI.Cf[block[i]].ndotl > 0)
	    sprintf(tmp,"%ld (%ld-%ld)",LI.Cf[block[i]].ndotl,
	        (long) firstdl[block[i]],
	        LI.Cf[block[i]].ndotl+firstdl[block[i]]-1);
	else
	    sprintf(tmp,"0");
	fprintf(f,"%-20s\n",tmp);
    }
    fprintf(f,"\n");

    /* Write mountains
       --------------- */
    if (opt_o & O_MOUNTAINS)
    {
	fprintf(f,"Mountains:\n");
	fprintf(f,"    No     Dim    Maximal Submountains\n");
	for (i = 0; i < bnmount; ++i)
	{
	    fprintf(f,"    %-7d%-7ld",i,bmdim[i]);
	    BsCopy(b,bsupof[i]);
	    BsClear(b,i);
	    for (k = 0; k < bnmount; ++k)
	    {
		if (!BsTest(b,k)) continue;
		BsMinus(b,bsupof[k]);
		BsSet(b,k);
	    }
	    for (k = 0; k < bnmount; ++k)
		if (BsTest(b,k)) fprintf(f,"%d ",k);
	    fprintf(f,"\n");
	}
	fprintf(f,"\n");
    }

    /* Write incidence matrix
       ---------------------- */
    if (opt_o & O_INCIDENCES)
    {
	MESSAGE(1,("  Incidence matrix (%d by %d)\n",bnmount,bnmount));
	fprintf(f,"Incidence matrix:\n");
	for (i = 0; i < bnmount; ++i)
	{
	    fprintf(f,"    %3d: ",i);
	    printbs(f,bsupof[i]);
	    fprintf(f,"\n");
	}
	fprintf(f,"\n");
    }

    /* Write dotted lines
       ------------------ */
    if (opt_o & O_DOTTEDLINES)
    {
	MESSAGE(1,("  Dotted lines (%d)\n",bndotl));
	fflush(stdout);
	fprintf(f,"Dotted lines:\n");
	for (i = 0; i < bndotl; ++i)
	{
		fprintf(f,"    ");
		printbs(f,bdotl[i]);
		fprintf(f,"\n");
	}
	fprintf(f,"\n");
    }

    /* Write submodules
       ---------------- */
    if (opt_o & O_SUBMODULES)
    {
	MESSAGE(1,("  Submodules (%d)\n",nsub));
	fflush(stdout);
    	g = openout(".sub");
	fprintf(f,"Submodules:\n");
	fprintf(f,"    No    Dim  Flags  Ident                           Max\n");
	for (i = 0; i < nsub; ++i)
	{
	    int *lp;

	    fprintf(f,"    %-6d%-5ld",i,subdim[i]);
	    putc(ismount[i] ? 'M' : ' ',f);
	    putc(israd[i] ? 'R' : ' ',f);
	    putc(issoc[i] ? 'S' : ' ',f);
	    fprintf(f,"    ");
	    k = printbs(f,sub[i]);
	    for (; k < 30; ++k) putc(' ',f);
	    fprintf(f,"  ");
	    for (lp = max[i]; *lp >= 0; lp += 2)
	    {
		if (lp != max[i]) putc(',',f);
		fprintf(f,"%d",*lp);
	    }
	    fprintf(f,"\n");
	    BsWrite(sub[i],g);
	}
	fprintf(f,"\n");
	fclose(g);
    }

    /* Radikal- und Sockelreihe
       ------------------------ */
    if (opt_o & O_RADICAL)
    {
	static int mult[LAT_MAXCF];
	long rdim;
	int layer;
	BitString_t
	    *rad = BsAlloc(bnmount),
	    *newrad = BsAlloc(bnmount),
	    *x = BsAlloc(bnmount),
	    *zero = BsAlloc(bnmount);

	MESSAGE(1,("  Radical series\n"));
	fprintf(f,"Radical series:\n");
	for (i = 0; i < bnmount; ++i) BsSet(rad,i);
	BsClearAll(zero);
	for (i = 0, rdim = 0; i < LI.NCf; ++i)
	    rdim += LI.Cf[i].dim * LI.Cf[i].mult;
	for (layer = 1; BsCompare(rad,zero) != 0; ++layer)
	{
	    BsClearAll(x);
	    BsClearAll(newrad);
	    MESSAGE(1,("Starting layer %d\n",layer));

	    /* Extend the zero module = x by all those mountains
	       which are contained in the radical and extend y
	       ------------------------------------------------- */
	    for (i = 0; i < bnmount && BsCompare(rad,x); ++i)
	    {
		if (BsTest(rad,i) && !(BsTest(x,i)))
		{
		    MESSAGE(2,("extend(%i)\n",i));
		    extend(x,i,0);
		    MESSAGE(2,("nextend(%i)\n",i));
		    extend(newrad,i,1);
		}
	    }

	    /* Find the irreducible factors in this layer
	       ------------------------------------------ */
	    memset(mult,0,sizeof(mult));
	    BsCopy(x,newrad);
	    for (i = 0; i < bnmount && BsCompare(rad,x); ++i)
	    {
	      if (BsTest(rad,i) && !(BsTest(x,i)))
	      {
		  extend(x,i,0);
		  k = isotype(i);
		  ++mult[k];
		  rdim -= LI.Cf[k].dim;
	      }
	    }
	    fprintf(f,"    Layer %d: Dim=%-4ld  ",layer,rdim);
	    for (i = 0; i < LI.NCf; ++i)
		for (; mult[i] > 0; --mult[i])
		    fprintf(f,"%s ",Lat_CfName(&LI,i));
	    fprintf(f,"\n");
	    BsCopy(rad,newrad);
	}
	fprintf(f,"\n");
    }
    fclose(f);


    if ((opt_o & O_EXTFILES) && (opt_o & O_SUBMODULES))
    {

	/* Write the .lat file
	   ------------------- */
	f = openout(".lat");
#if 0
	fprintf(f,"MeatAxe.IncidenceMatrix := [\n");
	for (i = 0; i < nsub; ++i)
	{
	    fprintf(f,"[");
	    for (k = 0; k < nsub; ++k)
	    {
		fprintf(f,"%d",BsIsSub(sub[i],sub[k]) ? 1 : 0);
		if (k < nsub-1)
		    fprintf(f,",");
	    }
	    fprintf(f,"]");
	    if (i < nsub-1)
		fprintf(f,",");
	    fprintf(f,"\n");
	}
	fprintf(f,"];\n");
#endif

	fprintf(f,"MeatAxe.Lattice := [\n");
	for (i = 0; i < nsub; ++i)
	{
	    int *lp = max[i];

	    fprintf(f,"[%ld,[",subdim[i]);
	    for (k = 0, lp = max[i]; *lp >= 0; lp += 2, ++k)
	    {
	    	fprintf(f,"[%d,%d]",lp[0]+1,lp[1]+1);
		if (lp[2] >= 0)
		{
		    fprintf(f,",");
		    if (k % 10 == 9) fprintf(f,"\n");
		}
	    }
	    if (i < nsub-1)
		fprintf(f,"]],\n");
	    else
		fprintf(f,"]]\n");
	}
	fprintf(f,"];\n");


	fclose(f);

	/* Write the .gra file
	   ------------------- */
	f = openout(".gra");
	fprintf(f,"%d\n",nsub);
	for (i = 0; i < nsub; ++i)
	{
	int *lp;

	putc(ismount[i] ? 'm' : '.',f);
	putc(israd[i] ? 'r' : '.',f);
	putc(issoc[i] ? 's' : '.',f);
	for (k = 0, lp = max[i]; *lp >= 0; lp += 2, ++k);
	fprintf(f," %2d",k);
	for (lp = max[i]; *lp >= 0; lp += 2)
	    fprintf(f," %d %d",lp[0],lp[1]);
	fprintf(f,"\n");
	}
	fclose(f);
    }
}



/* -----------------------------------------------------------------
   sameblock() - Test if irred. #i and #k belong to the same block
	(Actually, we check only if there is any mountain of
	irred. i contained in any mountain of irred. k and vice
	versa)
   ----------------------------------------------------------------- */

static int sameblock(int i, int k)

{
    int ii, kk;

    for (ii = firstm[i]; ii < firstm[i+1]; ++ii)
	for (kk = firstm[k]; kk < firstm[k+1]; ++kk)
	{
	    if (BsTest(xsubof[ii],kk)) return 1;
	    if (BsTest(xsubof[kk],ii)) return 1;
	}
    return 0;
}




/* -----------------------------------------------------------------
   nextblock() - Build next block; returns 0 if no more
   	constituents remain
   ----------------------------------------------------------------- */

static int nextblock()

{
    int i, k;

    ++blnum;
    for (i = 0; i < LI.NCf && done[i]; ++i);
    if (i >= LI.NCf) return 0;

    /* If -b was not used, build one single
       block containing all irreducibles.
       ------------------------------------ */
    if (!opt_b)
    {
	blsize = LI.NCf;
	for (i = 0; i < LI.NCf; ++i)
	{
	    block[i] = i;
	    done[i] = 1;
	}
	return 1;
    }

    /* Otherwise, make the next block
       ------------------------------ */
    MESSAGE(2,("Making next block (%d)\n",blnum));
    done[i] = 1;
    blsize = 1;
    block[0] = i;
    i = 0;
    while (i < blsize)
    {
	for (k = 0; k < LI.NCf; ++k)
    	    if (!done[k] && sameblock(block[i],k))
    	    {
	        done[k] = 1;
	        block[blsize++] = k;
	    }
	++i;
    }

    /* Sort block
       ---------- */
    MESSAGE(2,("Sorting\n"));
    for (i = 0; i < blsize; ++i)
	for (k = i+1; k < blsize; ++k)
	{
	    if (block[i] > block[k])
	    {
		int tmp = block[i];
		block[i] = block[k];
		block[k] = tmp;
	    }
	}
    if (MSG0)
    {
        printf("\nBlock %d: ",blnum);
	for (i = 0; i < blsize; ++i)
	    printf(" %s%s",LI.BaseName,Lat_CfName(&LI,block[i]));
	printf("\n");
        fflush(stdout);
    }
    return 1;
}



/* -----------------------------------------------------------------
   sort()
   ----------------------------------------------------------------- */

void sort()

{
    int i, k;
    BitString_t *x;

    MESSAGE(0,("Sorting\n"));
    for (i = 0; i < nsub; ++i)
    {
	for (k = i+1; k < nsub; ++k)
	{
	    if (BsIsSub(sub[k],sub[i]))
	    {
		x = sub[i];
		sub[i] = sub[k];
		sub[k] = x;
	    }
	}
    }
}

/* -----------------------------------------------------------------
   addtolist() - Find out if a given submodule is new. If yes, add
	the new sub module to the list
   ----------------------------------------------------------------- */

static void addtolist(BitString_t *x)

{
    int i;

    ++nadd;
    if (nadd % 1000 == 0)
    {
	MESSAGE(1,("Generation %d: %d tries, %d new submodules\n",
	    generation,nadd,nsub-oldnsub));
    }

    for (i = nsub-1; i >= 0; --i)
    {
	if (!BsCompare(sub[i],x))
	    return;
    }
    if (nsub >= MAXNSUB)
    {
	sort();
	findrsm();
	writeresult();	/* Write out what we have found so far */
	MTX_ERROR1("Too many submodules (> %d)",MAXNSUB);
    }
    sub[nsub] = BsAlloc(bnmount);
    BsCopy(sub[nsub++],x);
}





/* -----------------------------------------------------------------
   initblock()
   ----------------------------------------------------------------- */

static void initblock()

{
    int i, k, ii, kk, row, col;

    /* Find out the number of mountains in this block
       ---------------------------------------------- */
    bnmount = 0;
    for (i = 0; i < blsize; ++i)
	bnmount += LI.Cf[block[i]].nmount;

    /* Build the incidence matrix
       -------------------------- */
    MESSAGE(0,("Building incidence matrix\n"));
    fflush(stdout);
    for (i = 0; i < bnmount; ++i)
    {
	bsubof[i] = BsAlloc(bnmount);
	bsupof[i] = BsAlloc(bnmount);
    }
    row = 0;
    for (i = 0; i < blsize; ++i)
    {
        for (ii = firstm[block[i]]; ii < firstm[block[i]+1]; ++ii)
	{
	    col = 0;
	    for (k = 0; k < blsize; ++k)
	    {
		for (kk=firstm[block[k]]; kk<firstm[block[k]+1]; ++kk)
		{
		    if (BsTest(xsubof[ii],kk))
		    {
			BsSet(bsubof[row],col);
		        BsSet(bsupof[col],row);
		    }
		    ++col;
		}
	    }
	    bmdim[row] = xmdim[ii];
	    ++row;
	}
    }

    /* Build the dotted lines for one block
       ------------------------------------ */
    MESSAGE(0,("Building dotted lines\n"));
    fflush(stdout);
    bndotl = 0;
    for (i = 0; i < blsize; ++i)
    {
        for (ii = firstdl[block[i]]; ii < firstdl[block[i]+1]; ++ii)
	{
	    bdotl[bndotl] = BsAlloc(bnmount);
	    bdlspan[bndotl] = BsAlloc(bnmount);
	    col = 0;
	    for (k = 0; k < blsize; ++k)
	    {
		for (kk=firstm[block[k]]; kk<firstm[block[k]+1]; ++kk)
		{
		    if (BsTest(xdotl[ii],kk))
		    {
			BsOr(bdlspan[bndotl],bsupof[col]);
			BsSet(bdotl[bndotl],col);
		    }
		    ++col;
		}
	    }
	    ++bndotl;
	}
    }

    /* Initialize global variables
       --------------------------- */
    generation = 0;
    nsub = 0;		/* Number of submodules */
    lastgen = 0;	/* First submodule of previous generation */
    nadd = 0;
    oldnsub = 0;
    addtolist(BsAlloc(bnmount));	/* Null module */
}



/* -----------------------------------------------------------------
   cleanupblock() - Clean up after each block
   ----------------------------------------------------------------- */

static void cleanupblock()

{
    int i;

    for (i = 0; i < bnmount; ++i)
    {
	free(bsubof[i]);
	free(bsupof[i]);
    }
    for (i = 0; i < bndotl; ++i)
    {
	free(bdotl[i]);
    }
    if (opt_o & O_SUBMODULES)
    {
        for (i = 0; i < nsub; ++i)
        {
	    free(sub[i]);
	    free(max[i]);
        }
        free(ismount);
        free(israd);
        free(issoc);
        free(subdim);
        free(max);
    }
}



/* -----------------------------------------------------------------
   nextgen() - Make next generation (modules generated by n+1
	mountains)
   ----------------------------------------------------------------- */

void nextgen()

{
    int i, k, oldnsub = nsub;
    BitString_t *x;

    x = BsAlloc(bnmount);
    for (i = lastgen; i < oldnsub; ++i)
    {
	for (k = 0; k < bnmount; ++k)
	{
	    if (BsTest(sub[i],k))
		continue;
	    BsCopy(x,sub[i]);
	    extend(x,k,0);
	    addtolist(x);
	}
    }
    lastgen = oldnsub;
    ++generation;
}


static int SetFormatFlags(const char *c, int set)

{
    static int FirstTime = 1;
    if (FirstTime)
    {
	FirstTime = 0;
	opt_o = set ? 0 : O_ALL;
    }
    for (; *c != 0; ++c)
    {
	int flag = 0;
	switch (*c)
	{
	    case 'm': flag = O_MOUNTAINS; break;
	    case 's': flag = O_SUBMODULES; break;
	    case 'd': flag = O_DOTTEDLINES; break;
	    case 'e': flag = O_EXTFILES; break;
	    case 'r': flag = O_RADICAL; break;
	    case 'o': flag = O_SOCLE; break;
	    case 'i': flag = O_INCIDENCES; break;
	    default: 
		MTX_ERROR1("Unknown format flag '%c'",*c);
		return -1;
	}
	if (set)
	    opt_o |= flag;
	else
	    opt_o &= ~flag;
    }
    return 0;
}



static int ParseCommandLine()

{
    const char *c;

    /* Parse command line
       ------------------ */
    opt_G = AppGetOption(App,"-G --gap");
    if (opt_G) 
	MtxMessageLevel = -100;
    opt_b = AppGetOption(App,"-b --blocks");
    if ((c = AppGetTextOption(App,"-o --output",NULL)) != NULL)
    {
	if (SetFormatFlags(c,1) != 0)
	    return -1;
    }
    if ((c = AppGetTextOption(App,"-n --no-output",NULL)) != NULL)
    {
	if (SetFormatFlags(c,0) != 0)
	    return -1;
    }
    if (AppGetArguments(App,1,1) != 1)
	return -1;
    return 0;
}



static int Init(int argc, const char **argv)

{
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (ParseCommandLine() != 0)
	return -1;
    MESSAGE(0,("*** CALCULATE ALL SUBMODULES ***\n\n"));
    init(App->ArgV[0]);
    init2();
    return 0;
}





/* -----------------------------------------------------------------
   main()
   ----------------------------------------------------------------- */

int main(int argc, const char **argv)

{
    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return -1;
    }

    while (nextblock())
    {
	initblock();

	if (opt_o & O_SUBMODULES)
	{
	    do
	    {
		MESSAGE(0,("Generation %d: %d tr%s, %d new submodule%s\n",
		    generation,nadd,nadd == 1 ? "y":"ies",
		    nsub-oldnsub,nsub-oldnsub == 1 ? "" : "s"));
/*if (nsub > 100) exit(0);
*/
		nadd = 0;
		oldnsub = nsub;
		nextgen();
	    }
	    while (oldnsub != nsub);
	    sort();
	    findrsm();
	}
	else
	    MESSAGE(0,("Submodules not calculated\n"));

	writeresult();
	MESSAGE(0,("\n"));
	cleanupblock();
    }
    MtxCleanupLibrary();
    return 0;
}



/**
@page prog_mksub mksub - Find Submodules

<<<<<<< HEAD
@section syntax Command Line
=======
@section mksub_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
<pre>
mksub @em Options [-Gb] [-o @em Format] [-n @em Format] @em Name
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -G
  Produce output in GAP format.
@par -b
  Find blocks (see description).
@par -o @em Format
  Include the selected parts in the output (see description).
@par -n @em Format
  Exclude the selected parts from the output (see description).
@par @em Name
  Name of the representation.

<<<<<<< HEAD
@section inp Input Files
=======
@section mksub_inp Input Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
@par @em Name.cfinfo
  Constituent info file.
@par @em Name.inc
  Incidence matrix generated by @ref prog_mkinc "mkinc".
@par @em Name.dot
  Dotted-lines generated by @ref prog_mkdotl "mkdotl".
@par @em Name.mnt
  Mountain dimensions.

<<<<<<< HEAD
@section out Output Files
=======
@section mksub_out Output Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
@par @em Name.cfinfo
  Constituent info file.
@par @em Name.out
  Submodule lattice.
@par @em Name.lat
  Incidence matrix of the submodules (GAP).
@par @em Name.gra
  Submodule lattice for @ref prog_mkgraph "mkgraph".

<<<<<<< HEAD
@section desc Description
=======
@section mksub_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c

This program calculates the submodule lattice using the output
generated by @ref prog_mkinc "mkinc" and @ref prog_mkdotl "mkdotl".
In this final step no matrix operations are involved.
Instead the program works with bit strings representing the incidences
and dotted-lines.
The lattice may be decomposed into blocks using the "-b" option (see below).

Submodules are calculated generation by generation, the first
generation consisting of all submodules generated by one 
local submodule. In the n-th generation all submodules generated
by a submodule of generation n plus one local submodule are
calculated. If no more submodules appear, the algorithm terminates.

Output is written to three text files and one binary file. The first
file, @em Name.out, contains a list of irreducible constituents, 
the incidence matrix, the dimensions of all local submodules, a list of
dotted lines, a list of all submodules, the radical and socle series,
and a list of all mountains, i.e., local submodules.
The second output file is @em Name.lat. It contains the lattice as 
a list in GAP format. This list contains, for each
submodule, its dimension, maximal submodules and isomorphism types
of simple factors are given (see the example below). The third output 
file, @em Name.gra, contains a description
of the submodule lattice together with some additional information. 
This file is read by the @ref prog_mkgraph "mkgraph" program to produce a graphical 
representation of the lattice. A further output file, @em Name.sub 
contains the submodule lattice in binary format. This file is read by
the @ref prog_genmod "genmod" program.

@subsection blk Blocks
If the "-b" option is used, @b mksub tries to decompose the lattice
into blocks. By definition, a block is a set of one or more
composition factors which is closed under incidences of local
submodules. For a decomposition into blocks to be possible, there 
must be direct summands with no common irreducible constituent.
If a decomposition exists, the whole lattice can be reconstructed
from its blocks by forming direct sums.

Output files are created separately for each block, and a number
is appended to the name. For example, if the representation is
called "psl277" and the lattice decomposes into 3 blocks,
@b mksub creates 9 output files:
<pre>
psl211.out.1 psl211.lat.1 psl211.gra.1
psl211.out.2 psl211.lat.2 psl211.gra.2
psl211.out.3 psl211.lat.3 psl211.gra.3
</pre>


@subsection of Changing the Output Format

The default output as shown above may be changed by using the
"-o" option. @em Format is any combination of "m", "d", "i", "e", "s", "r", "s".
Each letter corresponds to a certain piece of output:
@par m
  Mountains.
@par d
  Dotted lines.
@par i
  Incidence matrix.
@par e
  ".lat" and ".gra" files.
@par s
  Submodules.
@par r
  Radical series.
@par o
  Socle series.

**/

