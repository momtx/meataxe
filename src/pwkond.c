////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Find peak words and condense
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <string.h>
#include <stdlib.h>




/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */
MTX_DEFINE_FILE_INFO

#define MAX_MODULES 50
#define MAXCF (3 * LAT_MAXCF)	    /* Max. number of different constituents */

int NumMods = 0;		    /* Number of modules */
static struct
{
    Lat_Info Info;		    /* Data from .cfinfo file */
    MatRep_t *Rep;		    /* Generators */
    WgData_t *Wg;		    /* Word generators */
    Matrix_t *SsBasis;		    /* Semisimplicity basis */
} ModList[MAX_MODULES];		    /* List of all modules */


int NumCf = 0;			    /* Number of inequivalent constituents */

/** Table of constituents
 **/
struct cf_struct
{
    MatRep_t *Gen;		    /* Generators */
    CfInfo *Info;		    /* Constituent info */
    WgData_t *Wg;		    /* Word generators */
    Lat_Info *Mod;		    /* Which module does this come from */
    int CfNo;			    /* Constituent index within that module */
    Matrix_t *PWNullSpace;	    /* Peak word null space */
    int Mult;			    /* In how many modules does it appear? */
    int CfMap[MAX_MODULES][2];	    /* (module,cf) */
} CfList[MAXCF];		    /* List of all constituents */



static int opt_G = 0;
static int opt_n = 0;			/* No condensation, PW only */
static int opt_p = 0;			/* Use full polynomials in PW search */
static int opt_t = 0;			/* Transform generators to std basis  */
static int opt_b = 0;			/* Calculate a semisimplicity basis */
static int opt_k = 0;                 /* Make kernel of PW */


#define MAXLOCK 100
static long include[MAXLOCK][2];
static int ninclude=0;
static long exclude[MAXLOCK][2];
static int nexclude=0;
static int PeakWordsMissing;		/* Number of missing peak words */


static MtxApplicationInfo_t AppInfo = { 
"pwkond", "Peakword Condensation",
"\n"
"SYNTAX\n"
"    pwkond [<Options>] <Name> [<Name> ...]\n"
"\n"
"ARGUMENTS\n"
"    <Name> .................. Name of the representation\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... GAP output (implies -Q)\n"
"    -n ...................... Find peak words only, do not condense\n"
"    -p ...................... Use full polynomials in peak word search\n"
"    -i <List> ............... Words to try first. Example: -i 100,20-35.\n"
"    -e <List> ............... Exclude words from search. Example: -e 3,20-99.\n"
"    -t ...................... Transform generators into standard basis\n"
"    -b ...................... Calculate a semisimplicity basis\n"
"    -k ...................... Compute kernel of peak words\n"

"\n"
"FILES\n"
"    <Name>.cfinfo ........... IO Constituent info file\n"
"    <Name>.{1,2,...} ........ I  Generators\n"
"    <Name><Cf>.{1,2...} ..... I  Generators on the constituents\n"
"    <Name><Cf>.{1,2...}k .... O  Condensed generators\n"
"    <Name><Cf>.{1,2...}.std   O  Condensed generators in std basis (with -t)\n"
"    <Name><Cf>.op ........... O  Spin-up script for standard basis (with -t)\n"
"    <Name><Cf>.np ........... O  Condensed peak word\n"
"    <Name><Cf>.im ........... O  Image used for condensation\n"
"    <Name><Cf>.k ............ O  Peakword kernel (with -k or without -n)\n"
"    <Name>.ssb .............. O  Semisimplicity basis (with -b)\n"
};



static MtxApplication_t *App = NULL;




/* -------------------------------------------------------------------------- 
   AddConstituent() - Add a constituent

   This function checks if a given constituent is already in <CfList>. If
   not, it is added to the list. If it is already in the list, <cf> is
   deleted.
   -------------------------------------------------------------------------- */

static int AddConstituent(MatRep_t *cf, CfInfo *info, int modno, int cfno)

{
    int i, m;
    for (i = 0; i < NumCf; ++i)
    {
	int rc;
	rc = IsIsomorphic(CfList[i].Gen,CfList[i].Info,cf,NULL,0);
	if (rc < 0)
	    return -1;
	if (rc) 
	    break;
    }

    if (i < NumCf)  /* Constituent was already in the list */
    {
	int m = CfList[i].Mult;
	MrFree(cf);
	CfList[i].CfMap[m][0] = modno;
    }
    else	    /* It's a new constituent */
    {
	CfList[i].Gen = cf;
	CfList[i].Info = info;
	CfList[i].Wg = WgAlloc(cf);
	CfList[i].Mod = &ModList[modno].Info;
	CfList[i].Mult = 0;
	++NumCf;
    }
    m = CfList[i].Mult;
    CfList[i].CfMap[m][0] = modno;
    CfList[i].CfMap[m][1] = cfno;
    CfList[i].Mult++;
    MESSAGE(1,("%s%s is constituent %d\n",ModList[modno].Info.BaseName,
	Lat_CfName(&ModList[modno].Info,cfno),i));
    return i;
}


/* -------------------------------------------------------------------------- 
   AddConstituents() - Add all constituents of a module

   This function calls <AddConstituent()> for each constituent of the 
   <mod>-th module in <ModList>. I.e., it adds all constituents of that
   module to the global constituent list and sets up the constituent map.
   -------------------------------------------------------------------------- */

static int AddConstituents(int mod)

{
    Lat_Info *li = &ModList[mod].Info;
    int i;
    for (i = 0; i < li->NCf; ++i)
    {
	int index;
	char fn[200];
	MatRep_t *cf;
	sprintf(fn,"%s%s",li->BaseName,Lat_CfName(li,i));
	cf = MrLoad(fn,li->NGen);
	if (cf == NULL)
	    return -1;
	index = AddConstituent(cf,li->Cf + i,mod,i);
	if (index < 0)
	    return -1;		    /* Error */
    }
    return 0;
}


/* --------------------------------------------------------------------------
   LoadConstituents() - Load the generators on all constituents.
   -------------------------------------------------------------------------- */

static int LoadConstituents()

{
    int i;

    for (i = 0; i < NumMods; ++i)
    {
	if (AddConstituents(i) < 0)
	{
	    MTX_ERROR1("Error while loading constituents for %s",
		ModList[i].Info.BaseName);
	    return -1;
	}
    }

    /* Sort the constituents by dimension to speed up the peak word search.
       -------------------------------------------------------------------- */
    for (i = 0; i < NumCf; ++i)
    {
	int k;
	for (k = i + 1; k < NumCf; ++k)
	{
	    if (CfList[i].Info->dim > CfList[k].Info->dim)
	    {
		struct cf_struct tmp = CfList[i];
		CfList[i] = CfList[k];
		CfList[k] = tmp;
	    }
	}
    }
    return 0;
}






static int IsCompatible(int i)

{
    if (   ModList[i].Info.NGen != ModList[0].Info.NGen
	|| ModList[i].Info.Field != ModList[0].Info.Field)
    {
	MTX_ERROR3("%s and %s: %d",App->ArgV[0],App->ArgV[i],
	    MTX_ERR_INCOMPAT);
	return 0;
    }
    return 1;
}



/* --------------------------------------------------------------------------
   LoadModules() - Load the .cfinfo files and generators of the modules.
   -------------------------------------------------------------------------- */

static int LoadModules()

{
    int i;

    /* Set the number of modules.
       -------------------------- */
    NumMods = App->ArgC;
    if (NumMods > MAX_MODULES)
    {
	MTX_ERROR1("Too many modules (max. %d allowed)",MAX_MODULES);
	return -1;
    }

    /* Read the .cfinfo files and load the generators (if needed).
       ----------------------------------------------------------- */
    for (i = 0; i < NumMods; ++i)
    {
	int k;
        if (Lat_ReadInfo(&ModList[i].Info,App->ArgV[i]) != 0)
	    return -1;
	if (!IsCompatible(i))
	    return -1;

	/* Clear any existing peak words.
	   ------------------------------ */
	for (k = 0; k < ModList[i].Info.NCf; ++k)
	    ModList[i].Info.Cf[k].peakword = -1;

	/* Read the generators, set up ss bases and word generators.
	   --------------------------------------------------------- */
	if (!opt_n || opt_k || opt_b)
	{
	    ModList[i].Rep = MrLoad(App->ArgV[i],ModList[i].Info.NGen);
	    ModList[i].Wg = WgAlloc(ModList[i].Rep);
	    if (opt_b)
	    {
		int dim = ModList[i].Rep->Gen[0]->Nor;
		ModList[i].SsBasis = MatAlloc(FfOrder,dim,dim);
	    }
	}
    }

    return 0;
}






/* ------------------------------------------------------------------
   gkond() - Generalized condensation of one matrix
   ------------------------------------------------------------------ */

static void gkond(const Lat_Info *li, int i, Matrix_t *b, Matrix_t *k, 
    Matrix_t *w, const char *name)

{
    char fn[LAT_MAXBASENAME+10];
    Matrix_t *x1, *x2;

    x1 = MatDup(k);
    MatMul(x1,w);
    x2 = QProjection(b,x1);
    sprintf(fn,"%s%s.%s",li->BaseName,Lat_CfName(li,i),name);
    MatSave(x2,fn);
    MatFree(x1);
    MatFree(x2);
}



/* ------------------------------------------------------------------
   Standardize() - Bringe Kompositionsfaktor auf Standardform.
   Schreibt die neuen Erzeuger in XXX.std.n und die Operationen in
   XXX.op.
   ------------------------------------------------------------------ */

static void Standardize(int cf)

{
    int k, m;
    Matrix_t *sb;
    IntMatrix_t *script = NULL;
    Matrix_t *std[MAXGEN];

    /* Make the spin-up script for the standard basis and
       transform the generators.
       -------------------------------------------------- */
    MESSAGE(0,("  Transforming to standard basis\n"));
    sb = SpinUp(CfList[cf].PWNullSpace,CfList[cf].Gen,
	SF_FIRST|SF_CYCLIC|SF_STD,&script,NULL);
    ChangeBasisOLD(sb,CfList[cf].Gen->NGen,
	(const Matrix_t **)CfList[cf].Gen->Gen,std);
    MatFree(sb);

    /* Write the transformed generators and the spin-up script.
       -------------------------------------------------------- */
    for (m = 0; m < CfList[cf].Mult; ++m)
    {
	char fn[200];
	Lat_Info *li = &ModList[CfList[cf].CfMap[m][0]].Info;
	int i = CfList[cf].CfMap[m][1];
	sprintf(fn,"%s%s.op",li->BaseName,Lat_CfName(li,i));
	MESSAGE(2,("Write operations to %s\n",fn));
	if (ImatSave(script,fn) != 0)
	    MTX_ERROR("Cannot write .op file");
	for (k = 0; k < li->NGen; ++k)
	{
    	    sprintf(fn,"%s%s.std.%d",li->BaseName,Lat_CfName(li,i),k+1);
    	    MESSAGE(2,(" %s",fn));
	    MatSave(std[k],fn);
	}
    }

    for (k = 0; k < ModList[0].Info.NGen; ++k)
	MatFree(std[k]);
    ImatFree(script);
}



/* ------------------------------------------------------------------
   CfPosition() - Find the starting point for a constituent.
   ------------------------------------------------------------------ */

static int CfPosition(const Lat_Info *li, int cf)

{
    int pos = 0;
    int i;
    MTX_VERIFY(cf >= 0 && cf < li->NCf);
    for (i = 0; i < cf; ++i)
	pos += li->Cf[i].dim * li->Cf[i].mult;
    return pos;
}


/* ------------------------------------------------------------------
   kond() - Generalized condensation for one irreducible
   ------------------------------------------------------------------ */

static void kond(int mod, int cf)

{
    const Lat_Info *li = &ModList[mod].Info;
    char fn[LAT_MAXBASENAME+10];
    Matrix_t *peakword, *kern, *m, *k, *pw;
    int j, pwr;
		
    /* Make the peak word, find its stable power, 
       and calculate both kernel and image.
       ------------------------------------------ */
    peakword = WgMakeWord(ModList[mod].Wg,li->Cf[cf].peakword);
    MatInsert_(peakword,li->Cf[cf].peakpol);
    pw = MatDup(peakword);
    StablePower_(peakword,&pwr,&kern);
    MESSAGE(0,("pwr=%d, nul=%d, ",pwr,kern->Nor));
    if (kern->Nor != li->Cf[cf].mult * li->Cf[cf].spl)
	MTX_ERROR("Something is wrong here!");
    MatEchelonize(peakword);

    /* Write out the image
       ------------------- */
     if (!opt_n)
     {
	sprintf(fn,"%s%s.im",li->BaseName,Lat_CfName(li,cf));
	MatSave(peakword,fn);
     }

    /* Write out the `uncondense matrix'
       --------------------------------- */
    m = QProjection(peakword,kern);
    k = MatInverse(m);
    MatFree(m);
    MatMul(k,kern);
    sprintf(fn,"%s%s.k",li->BaseName,Lat_CfName(li,cf));
    MatSave(k,fn);

    /* Condense all generators
       ----------------------- */
    MESSAGE(1,("("));
    for (j = 0; j < li->NGen; ++j)
    {
	sprintf(fn,"%dk",j+1);
	gkond(li,cf,peakword,k,ModList[mod].Rep->Gen[j],fn);
	MESSAGE(1,("%d",j+1));
    }
    MESSAGE(1,(")"));

    /* Condense the peak word
       ---------------------- */
    gkond(li,cf,peakword,k,pw,"np");

    /* Calculate the semisimplicity basis.
       -----------------------------------  */
    if (opt_b)
    {
	Matrix_t *seed, *partbas;
	int pos = CfPosition(li,cf);
	seed = MatNullSpace_(pw,0);
	partbas = SpinUp(seed,ModList[mod].Rep,SF_EACH|SF_COMBINE|SF_STD,NULL,NULL);
        MatFree(seed);
	MESSAGE(0,(", %d basis vectors",partbas->Nor));
	if (MatCopyRegion(ModList[mod].SsBasis,pos,0,partbas,0,0,-1,-1) != 0)
	{
	    MTX_ERROR1("Error making basis - '%s' is possibly not semisimple",
		li->BaseName);
	}
        MatFree(partbas);
    }
    MatFree(pw);

    MESSAGE(0,("\n"));

    MatFree(k);
    MatFree(kern);
    MatFree(peakword);

}




static void Condense(int cf)
{
    int k;
    for (k = 0; k < CfList[cf].Mult; ++k)
    {
	int mod = CfList[cf].CfMap[k][0];	
	int i = CfList[cf].CfMap[k][1];	

	MESSAGE(0,("  Condensing %s%s: ",ModList[mod].Info.BaseName,
	    Lat_CfName(&ModList[mod].Info,i)));
	kond(mod,i);
    }
}


static String GapPrintPoly(const Poly_t *pol)
{ 
    String s = StrAlloc(10);
    int i;
    StrAppend(&s, "[");
    for (i = 0; i < pol->Degree; ++i)
	StrAppendF(&s,"%s,",FfToGap(pol->Data[i]));
    StrAppendF(&s,"%s]",FfToGap(pol->Data[i]));
    return s;
}


static String GapPrintWord(const WgData_t *b, long n)
{
    String buffer = StrAlloc(10);

    WgDescribeWord((WgData_t *)b, n);
    int *x;
    StrAppend(&buffer, "[");
    for (x = b->Description; *x != -1; ) {
	int first = 1;
	if (x != b->Description) StrAppend(&buffer, ",");
	StrAppend(&buffer, "[");
	do {
	    long gen = *x++;
	    if (!first) StrAppend(&buffer, ",");
	    first = 0;
	    StrAppendF(&buffer,"%d", gen + 1);
	} while (*x != -1);
	StrAppend(&buffer, "]");
	++x;
    }
    StrAppend(&buffer, "]");
    return buffer;
}




/* --------------------------------------------------------------------------
   WriteOutput() - Write all output files

   Description:
     This function writes the .cfinfo file and the semisimplicity basis. It
     is called after a new peak word has been found.
   -------------------------------------------------------------------------- */

static void WriteOutput(int final)

{
    int i;
    for (i = 0; i < NumMods; ++i)
    {
	Lat_WriteInfo(&ModList[i].Info);
	if (opt_b)
	{
	    char fn[200];
	    sprintf(fn,"%s.ssb",ModList[i].Info.BaseName);
	    MESSAGE(1,("Writing semisimplicity basis to %s\n",fn));
	    MatSave(ModList[i].SsBasis,fn);
	}
    }
    if (!final)
	return;

    /* Write GAP output.
       ----------------- */
    if (opt_G)
    {
	int m, i;
	printf("MeatAxe.PeakWords := [\n");
	for (m = 0; m < NumMods; ++m)
	{
            const Lat_Info * const ModInfo = &ModList[m].Info;
	    printf("# module: %s\n[\n", ModInfo->BaseName);
	    for (i = 0; i < ModInfo->NCf; ++i)
	    {
    	        const CfInfo * const Cf = ModInfo->Cf + i;
		printf("    # irreducible factor: %s\n", Lat_CfName(ModInfo,i) );
		String ws = GapPrintWord(CfList[i].Wg,Cf->peakword);
		String ps = GapPrintPoly(Cf->peakpol);
		printf("    [ %ld, %s, %s ]%s\n", Cf->peakword, ws.S, ps.S,
		     i == ModInfo->NCf - 1 ? "" : "," );
		StrFree(&ws);
		StrFree(&ps);
	    }
	    printf(m < NumMods - 1 ? "],\n" : "]\n");
	}
	printf("];\n");
    }

}




/* --------------------------------------------------------------------------
   CopyPeakWordToAllModules() - Copy the peak word and polynomial just found 
   to all modules having an appropriate constituent. Also, print a nice 
   message showing the peak word and to which constituents it applies.
   -------------------------------------------------------------------------- */

static void CopyPeakWordToAllModules(int i)

{
    CfInfo *info = CfList[i].Info;
    const Poly_t *pp = info->peakpol;
    const int pw = info->peakword;
    int k;

    MESSAGE(0,("Peak word for"));
    for (k = 0; k < CfList[i].Mult; ++k)
    {
	int mod = CfList[i].CfMap[k][0];
	int l = CfList[i].CfMap[k][1];
	MESSAGE(0,("%c%s%s",k == 0 ? ' ': ',', 
	    ModList[mod].Info.BaseName,Lat_CfName(&ModList[mod].Info,l)));
	/* Copy peak word and peak polynomial to the other modules */
	if (k > 0)
	{
	    ModList[mod].Info.Cf[l].peakword = pw;
	    ModList[mod].Info.Cf[l].peakpol = PolDup(pp);
	}
    }
    MESSAGE(0,(" is %d (%s)",pw,WgSymbolicName(CfList[i].Wg,pw)));
    if (MSG0)
	PolPrint(", pol",pp);
}



/* --------------------------------------------------------------------------
   PeakWordFound() - A peak word has been found

   Arguments:
     <w>: The word number

   Description:
     This function is called each time a peak word is found. Depending on the
     command line options we condense the peak word, and transform the
     generators to standard basis.
     We also write the .cfinfo file each time. This allows the user to kill
     the running program and continue with the peak words found so far.
   -------------------------------------------------------------------------- */

static void PeakWordFound(int i)

{
    CopyPeakWordToAllModules(i);
    if (!opt_n || opt_k)
	Condense(i);		/* Condense */
    if (opt_t)
	Standardize(i);		/* Transform to std basis */
    WriteOutput(0);
    --PeakWordsMissing;
}




static int isexcluded(long w)

{	
    int i;

    for (i = 0; i < nexclude; ++i)
	if (w >= exclude[i][0] && w <= exclude[i][1])
	    return 1;
    return 0;
}


static void parselist(const char *c, long list[][2], int *count)

{	
    long a, b;

    while (*c != 0)
    {	
	a = b = 0;
	while (*c >= '0' && *c <= '9')
	    a = a * 10 + (*c++ - '0');
	if (*c == '-')
	{	
	    ++c;
	    while (*c >= '0' && *c <= '9')
		b = b * 10 + (*c++ - '0');
	}
	else
	    b = a;
	if (a == 0 || b == 0 || a > b) 
	    MTX_ERROR("BAD ARGUMENTS");
	list[*count][0] = a;
	list[*count][1] = b;
	++*count;
	if (*c == ',') 
	    ++c;
    }
}



static int ParseCommandLine()

{
    const char *c;

    opt_G = AppGetOption(App,"-G --gap");
    opt_n = AppGetOption(App,"-n --no-condensation");
    opt_p = AppGetOption(App,"-p --use-polynomials");
    opt_t = AppGetOption(App,"-t --make-std-basis");
    opt_b = AppGetOption(App,"-b --make-ss-basis");
    opt_k = AppGetOption(App,"-k --make-pw-kernel");
    while ((c = AppGetTextOption(App,"-e --exclude",NULL)) != NULL)
	parselist(c,exclude,&nexclude);
    while ((c = AppGetTextOption(App,"-i --include",NULL)) != NULL)
	parselist(c,include,&ninclude);
    if (AppGetArguments(App,1,MAX_MODULES) < 0)
	return -1;
    if (opt_G) 
	MtxMessageLevel = -100;
    return 0;
}






static int Init(int argc, const char **argv)

{
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (ParseCommandLine() != 0)
	return -1;
    MESSAGE(0,("*** PEAK WORD CONDENSATION ***\n\n"));

    if (LoadModules() != 0)
	return -1;
    if (LoadConstituents() != 0)
	return -1;
    PeakWordsMissing = NumCf;

    return 0;
}




/* ------------------------------------------------------------------
   try() - Try another word
   ------------------------------------------------------------------ */

static void addid(Matrix_t *m, FEL f)
{
    PTR x;
    long i;
    if (f == FF_ZERO)
	return;
    FfSetNoc(m->Noc);
    for (i = 0, x = m->Data; i < m->Nor; ++i, FfStepPtr(&x))
	FfInsert(x,i,FfAdd(FfExtract(x,i),f));
}


static int try2(long w, FEL f)

{
    int i, ppos = -1;
    long nul;

    MESSAGE(3,("Word %ld+%dI:",w,FfToInt(f)));
    for (i = 0; i < NumCf; ++i)  /* For each composition factor... */
    {
	Matrix_t *word;
	
	word = WgMakeWord(CfList[i].Wg,w);
	addid(word,f);
	nul = MatNullity__(MatDup(word));
        MESSAGE(3,(" %ld",nul));
	if (nul != 0 && nul != CfList[i].Info->spl)
	{
	    MatFree(word);
	    MESSAGE(3,("failed\n"));
	    return -1;
	}
	if (nul == CfList[i].Info->spl)
	{
	    if (ppos >= 0 || CfList[i].Info->peakword > 0)
	    {
		MatFree(word);
	    	MESSAGE(3,("failed\n"));
		return -1;  /* Nullity should be 0 */
	    }
	    nul = MatNullity__(MatMul(MatDup(word),word));
	    if (nul != CfList[i].Info->spl)
	    {
		MatFree(word);
	    	MESSAGE(3,("failed (nullity not stable)\n"));
		return -1;  /* Nullity is not stable */
	    }
	    ppos = i;
	}
        MatFree(word);
    }
    MESSAGE(3,("\n"));

    if (ppos > -1) /* we have found a new peak word */
    {
	Poly_t *pp;
	Matrix_t *word;
	CfList[ppos].Info->peakword = w;

	/* Rechne nochmals den Nullraum aus (wir brauchen
	   ihn spaeter fuer die Standardform.
	   ----------------------------------------------- */
	word = WgMakeWord(CfList[ppos].Wg,w);
	addid(word,f);
	CfList[ppos].PWNullSpace = MatNullSpace__(word);

	/* Erzeuge das Peakpolynom (hier immer vom Grad 1)
	   ----------------------------------------------- */
	pp = PolAlloc(FfOrder,1);
	pp->Data[0] = f;
	CfList[ppos].Info->peakpol = pp;
        PeakWordFound(ppos);
    }
    return ppos;
}





/* --------------------------------------------------------------------------
   Try_Linear() - Find peak word using linear polynomials

   Arguments:
     <w>: The word number

   Description:
     This function checks if there is one or more peak words $W_w+\lambda 1$
     with $\lambda\in F$.
   -------------------------------------------------------------------------- */

static void Try_Linear(long w)

{
    long f;
    for (f = 0; f < FfOrder && PeakWordsMissing > 0; ++f)
	try2(w,FfFromInt(f));
}



/* ------------------------------------------------------------------
   tryp2()
   ------------------------------------------------------------------ */

static int tryp2(long w, int cf, Poly_t *pol)

{
    int i;

    for (i = 0; i < NumCf; ++i)
    {
	Matrix_t *word, *wordp;
	long nul;

	if (i == cf)
	    continue;
	word = WgMakeWord(CfList[i].Wg,w);
	wordp = MatInsert(word,pol);
	MatFree(word);
	nul = MatNullity__(wordp);
	if (nul != 0) 
	    return -1;
    }
    return 0;
}





/* ------------------------------------------------------------------
   try_p() - Try another word using full minimal polynomials
   ------------------------------------------------------------------ */

static int try_p(long w)

{
    int i;

    for (i = 0; i < NumCf; ++i)
    {
	Matrix_t *word;
	FPoly_t *mp;
	int k;

	if (CfList[i].Info->peakword > 0)  
	    continue;			/* We already have a peak word */
	word = WgMakeWord(CfList[i].Wg,w);
	mp = MinPol(word);
        if (MSG3)
	{
	    printf("Constituent %d, minpol =\n",i);
	    FpPrint(NULL,mp);
        }
	for (k = 0; k < (int) mp->NFactors; ++k)
	{
	    if (mp->Factor[k]->Degree * mp->Mult[k] == CfList[i].Info->spl)
	    {
	       Matrix_t *wp, *wp2;
	       long nul;

	       if (MSG3)
	       {
		   printf("%d, ",i);
		   PolPrint("factor",mp->Factor[k]);
	       }
	       if (tryp2(w,i,mp->Factor[k]) == -1) 
		   continue;

	       /* Check if the nullity is stable
	          ------------------------------ */
	       wp = MatInsert(word,mp->Factor[k]);
	       wp2 = MatMul(MatDup(wp),wp);
	       MatFree(wp);
	       nul = MatNullity__(wp2);
	       if (nul != CfList[i].Info->spl) 
		   continue;
	       break;
	    }
	}

	if (k < (int) mp->NFactors)
	{
	    CfList[i].Info->peakword = w;
	    CfList[i].Info->peakpol = PolDup(mp->Factor[k]);
	    CfList[i].PWNullSpace = MatNullSpace__(MatInsert(word,mp->Factor[k]));
            PeakWordFound(i);
	}
	else
	    k = -1;	/* Not found */
	FpFree(mp);
	MatFree(word);
	if (k >= 0)
	    return i;
    }
    return -1;
}



static void TryNext(long w)

{
    static int count = 0;

    if (isexcluded(w))
	return;
    if ((MSG1 && (count % 50 == 0)) || MSG2)
    {
	printf("Word %ld\n",w);
	fflush(stdout);
    }
    if (opt_p) try_p(w); else Try_Linear(w);
}






/* --------------------------------------------------------------------------
   main()
   -------------------------------------------------------------------------- */

int main(int argc, const char **argv)

{
    int i;
    long w;

    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }


    /* Find peak words
       --------------- */
    for (i = 0; PeakWordsMissing > 0 && i < ninclude; ++i)
    {
	if (i == 0)
	    MESSAGE(1,("Trying words from inclusion list\n"));
	for (w = include[i][0]; w <= include[i][1]; ++w)
	    TryNext(w);
    }
    for (w = 1; PeakWordsMissing > 0; ++w)
	TryNext(w);


    WriteOutput(1);
    if (App != NULL) AppFree(App);
    return 0;
}


/**
@page prog_pwkond pwkond - Peak Word Condensation

@section pwkond_syntax Command Line
<pre>
pwkond @em Options [-Gnptbk] [-i @em List] [-e @em List] @em Name [@em Name ...]
</pre>

@par @em Options
Standard options, see @ref prog_stdopts
@par -G
Produce output in GAP format.
@par -n
Find peak words only, do not condense.
@par -p
Use full polynomials instead of linear factors only in peak word search.
@par -t
Transform generators into standard basis.
@par -b
Calculate a semisimplicity basis.
@par -k
Compute kernel of peak words.
@par -i @em List
Word to try first, for example "-i 100,20-35".
@par -e @em List
Words to be excluded, for example "-e 3,20-99".

@par @em Name
Name of the representation.

@section pwkond_inp Input Files
@par @em Name.cfinfo
Constituent info file.
@par @em Name.1, @em Name.2, ...
Generators.

@par @em NameCF.1, @em NameCF.2, ...
Generators on the irreducible constituents. NameCF is the representation name
followed by the dimension and a letter to distinguish inequivalent constituents
of the same dimension (see @ref prog_chop "chop").

@section pwkond_out Output Files
@par @em Name.cfinfo
Constituent info file.

@par @em NameCF.1k, @em NameCF.2k, ...
Condensed generators on the irreducible constituents.

@par @em NameCF.1.std, @em NameCF.2.std, ...
Condensed generators in standard basis (with -t).

@par @em NameCF.op
Spin-up script for the standard basis (with -t).

@par @em NameCF.np
Condensed peak word.

@par @em NameCF.im
Image used for the condensation.

@par @em NameCF.k
Peak word kernel (with -k or without -n).

@par @em NameCF.ssb
Semisimplicity basis (with -b).


@section pwkond_desc Description

After the irreducible constituents of a module, or a number of modules,
have been found with @ref prog_chop "chop", this program can be used
- to calulate peak words for the constituents,
- to condense the module using the peak words,
- to transform the generators on the constituents to the standard
  basis as defined by the peak word kernel, and
- to calculate a basis reflecting the direct decomposition of the 
  module, if the module is semisimple.
By definition, a "peak word" for the i-th constituent is an algebra element which has
minimal nullity on the i-th constituent and which operates regularly (i.e., with nullity 0)
on the other constituents.
Als for identifying words (see @ref prog_chop "chop"), the nullity of a peak word on
its constituent is equal to the degree of the splitting field for that
constituent.

When more than one module is specified on the command line, the peak words found by
@b pwkond are "global", i.e., each peak word selects 
exactly one of the constituents of alle the modules. Running @b pwkond 
successively on two modules does not generally produce global peak 
words, since a peak word found for module M may have a non-zero
nullity on a different constituent that occurs in another module N
but not in M.

The -e option can be used to exclude certain words from the search. 
@em List is a list of integers or ranges of integers, for example
"-e 57,82-112,289".
Using "-i" you can specify a list of words which will be tested first. 
This can significantly reduce computation time if you already know one 
or more peak words for a given module. 
The "-n" option disables the condensation phase. If this option is used,
the program stops after the peak words have been found.
If the "-t" option is specified, @b pwkond transforms the generators of all
irreducible constituents to the standard basis defined by the peak word.

For each composition factor there are several output files. If, for
example, one composition factor is X10a, @b pwkond will produce
the following files:
- X10a.std.1 and X10a.std.1 are the operation of the generators on
  the constituent with respect to the standard basis defined by the
  peak word. These files are created only if the `-t' option is used.
- X10a.op Spin-up script for the standard basis. See ZSB for details.
- X10a.1k and X10a.2k are th action of the generators on the 
  condensed module.
- X10a.np Condensed peak word. This is a nilpotent matrix.
- X10a.im Image of the peak word.
- X10a.k Kernel of the peak word.
The .cfinfo file is written each time a peak word is found. So, if
the program does not terminate or dies unexpectedly the information about
the peak words found so far are not lost.

If the module is semisimple, @b pwkond can 
calculate a basis that respects the decomposition into irreducible
constituents. With respect to this basis, the generators are in block
diagonal form, where the blocks occur in the order determined by @ref prog_chop "chop".
All blocks corresponding to the same constituent are equal, not only 
equivalent, and the blocks occur in their "natural" order (as defined by
@ref prog_chop "chop"). This is essential for the tensor condensation procedure 
(see @ref prog_precond "precond"). To calculate the semisimplicity basis, use the
"-b" option. 
The basis is written to @em Name.ssb. Using "-b" with a module that is not
semisimple produces undefined results. Most probably, @b pwkond will stop 
with the error message "row index out of range", or it will write a
singular matrix to @em Name.ssb.

@section pwkond_impl Implementation Details
Internally, a peak word is represented by a pair (n,p) where n is
the canonical number of the word (See @ref prog_zmw "zmw"), and p is a 
polynomial. The peak word represented by this pair is p(Wn), Wn 
being the n-th word. Without "-p", @b pwkond considers only linear 
polynomials. If the "-p" option is used, @b pwkond can find polynomials 
of any degree.

Whenever a peak word is found, the generalized condensation 
is calculated as follows: The peakword is caculated as a matrix acting on V,
which is then repeatedly raised to higher powers until the nullity stabilizes.
The stable nullity equals the multiplicity k of the constituent times the
degree [E:F] of the splitting field extension.
Having a power w^N of the peakword with stable nullity,
the condensation onto its kernel, i.e., the projection of V onto V/w^N(V),
is determined in the same way as in the @ref prog_zqt "zqt" program.
**/

