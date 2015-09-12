////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Chop a representation (Calculate composition series)
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <meataxe.h>
#include <string.h>
#include <stdlib.h>


/* --------------------------------------------------------------------------
   Some constants
   -------------------------------------------------------------------------- */

#define MAXTRIES 10000	    /* Number of tries before FindIdWord() fails */
#define MAXFP	6	    /* Fingerprint size */

MTX_DEFINE_FILE_INFO 



/* --------------------------------------------------------------------------
   typedefs
   -------------------------------------------------------------------------- */

/** A submodule. */
typedef struct nodestruct
{
    struct nodestruct *sub, *quot, *parent;
    int dim, num;		/* Dimension and number */
    MatRep_t *Rep;		/* Generators */
    MatRep_t *TrRep;		/* Transposed generators */
    long mult;			/* Multiplicity */
    long spl;			/* Degree of splitting field */
    long idword;		/* Word used for std basis */
    Poly_t *idpol;		/* Polynomial "  "    "    */
    Poly_t *f1, *f2;		/* characteristic Polynomial c=f1*f2*/
    int fprint[MAXFP];		/* Fingerprint */
    Set_t *badwords;
    Matrix_t *nsp, *nsptr;	/* Null space */
    Matrix_t *subsp;		/* Subspace */
    long ggt;			/* g.c.d. of nullities */
    WgData_t *wg;		/* Used by the word generator */
    long wnum;
    Matrix_t *word;
} node_t;


/* --------------------------------------------------------------------------
   Function prototypes
   -------------------------------------------------------------------------- */

static int Chop(node_t *n);


/* --------------------------------------------------------------------------
   Global variables
   -------------------------------------------------------------------------- */

node_t *root;		    /* Root node of the constituent tree */
long opt_deglimit = -1;	    /* Max. degree of irred. polynomials */
long opt_nullimit;	    /* Max nullity for exhaustive vector search */
node_t *irred[LAT_MAXCF];   /* List of irred. constituents */
long firstword = 1;
int opt_G = 0;		    /* GAP output */
int opt_i = 0;		    /* -i: read an existing .cfinfo file */
Set_t *goodwords;	    /* List of `good' words */
Lat_Info LI;		    /* Data for .cfinfo */

static long stat_svsplit = 0;	/* Statistics */
static long stat_cpirred = 0;
static long stat_nssplit = 0;
static long stat_dlsplit = 0;
static long stat_irred = 0;
static long stat_exsplit = 0;



static MtxApplicationInfo_t AppInfo = { 
"chop", "Find irreducible constituents",
"SYNTAX\n"
"    chop [<Options>] <Name>\n"
"\n"
"ARGUMENTS\n"
"    <Name> .................. Name of the representation\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... GAP output (implies -Q)\n"
"    -g <NGen> ............... Set number of generators (default is 2)\n"
"    -s <Word> ............... Start with word number <Word>\n"
"    -n <MaxNul> ............. Set limit on nullity\n"
"    -d <MaxDeg> ............. Set limit on degrees of polynomials\n"
"    -i ...................... Read <Name>.cfinfo, if it exists\n"
"\n"
"FILES\n"
"    <Name>.{1,2,...} ........ I Generators\n"
"    <Name>.cfinfo ........... O Constituent info file\n"
"    <Name><Cf>.{1,2...} ..... O Generators on the constituents\n"
};

static MtxApplication_t *App = NULL;



/* ------------------------------------------------------------------
   CreateNode() - Create a new node
   ------------------------------------------------------------------ */

static node_t *CreateNode(MatRep_t *rep, node_t *parent)

{
    node_t *n = ALLOC(node_t);
    memset(n,0,sizeof(node_t));

    n->parent = parent;
    n->Rep = rep;
    n->dim = rep->Gen[0]->Nor;
    n->num = -1;
    n->mult = -1;
    n->spl = -1;
    n->idword = -1;
    n->badwords = SetAlloc();
    n->nsp = n->nsptr = NULL;
    n->wg = WgAlloc(n->Rep);
    if (n->wg == NULL) 
    {
	MTX_ERROR("WgAlloc() failed");
	return NULL;
    }
    n->wnum = -1;
    return n;
}

#define MATFREE(x) { if ((x) != NULL) { MatFree(x); (x)=NULL; } }

/* --------------------------------------------------------------------------
   CleanUpNode() - Clean up a node_t structure

   Arguments:
     <n>: Pointer to the <node_t> structure
     <complete>: Free generators, too.

   Description:
     This function releases all memory allocated internally in a <node_t>
     structure. If <complete> is zero, the generators are not freed.
   -------------------------------------------------------------------------- */

static void CleanUpNode(node_t *n, int complete)

{
    if (complete && n->Rep != NULL)
    {
	MrFree(n->Rep);
	n->Rep = NULL;
    }
    if (n->TrRep != NULL)
    {
	MrFree(n->TrRep);
	n->TrRep = NULL;
    }
    MATFREE(n->subsp);
    if (n->wg != NULL) { WgFree(n->wg); n->wg = NULL; }
    MATFREE(n->nsp);
    MATFREE(n->nsptr);
    MATFREE(n->word);
    if (n->f1 != NULL) { PolFree(n->f1); n->f1 = NULL; }
    if (n->f2 != NULL) { PolFree(n->f2); n->f2 = NULL; }
}




/* ------------------------------------------------------------------
   MakeWord() - Make a word
   ------------------------------------------------------------------ */

static void MakeWord(node_t *n, long w)

{
    MATFREE(n->word);
    n->word = WgMakeWord(n->wg,w);
    n->wnum = w;
}



/* --------------------------------------------------------------------------
   InsertWord() - Insert the current word into a polynomial and calculate 
   the kernel of both the resulting matrix (nsp), and of the transposed
   matrix (nsptr).
   -------------------------------------------------------------------------- */

static void InsertWord(node_t *n, Poly_t *p)

{
    Matrix_t *m, *mt;

    m = MatInsert(n->word,p);
    mt = MatTransposed(m);
    MATFREE(n->nsp);
    n->nsp = MatNullSpace__(m);
    MATFREE(n->nsptr);
    n->nsptr = MatNullSpace__(mt);
}


/* --------------------------------------------------------------------------
   CreateRoot() - Read the generators and create the root node of the 
   constituent tree.
   -------------------------------------------------------------------------- */

static int CreateRoot()

{
    MatRep_t *rep;
    FILE *f;
    char fn[200];

    if (opt_i)
    {
	/* Try to get the number of generators from the cfinfo file.
	   --------------------------------------------------------- */
	sprintf(fn,"%s.cfinfo",LI.BaseName);
	f = SysFopen(fn,FM_NOERROR|FM_READ);
	if (f != NULL)
	{
	    Lat_Info info;
    	    fclose(f);
	    Lat_ReadInfo(&info,LI.BaseName);
	    LI.NGen = info.NGen;
	    if (info.NGen > MAXGEN)
		MTX_ERROR1("Too many generators (max. %d -- "
		    "configured in meataxe.h)",MAXGEN);
	    MESSAGE(1,("Set number of generators = %d from %s\n",info.NGen,fn));
	}
    }

    rep = MrLoad(LI.BaseName,LI.NGen);
    LI.Field = FfOrder;
    if (rep == NULL)
	return -1;
    root = CreateNode(rep,NULL);
    if (root == NULL)
	return -1;
    return 0;
}



/* ------------------------------------------------------------------
   IsBadWord() - Look up a word in the 'badwords' list
   ------------------------------------------------------------------ */

static int IsBadWord(long w,node_t *n)

{
    if (n == NULL) 
	return 0;
    if (n->badwords != NULL && SetContains(n->badwords,w))
	return 1;
    return IsBadWord(w,n->parent);
}



/* ------------------------------------------------------------------
   PrintCompositionSeries() - Write the composition series to stdout
   ------------------------------------------------------------------ */

static void PrintCompositionSeries(node_t *n)

{
    static int count = 0;	/* Characters written in one line */
    int i;

    if (n == NULL)
    {
	MTX_ERROR1("n=NULL: %E",MTX_ERR_BADARG);
	return;
    }
    if (n->sub == NULL)		/* Irreducible */
    {
	if (count >= 20)
	{   
	    printf("\n");
	    count = 0;
	}
	for (i = 0; LI.Cf[i].dim != n->dim || LI.Cf[i].num != n->num; ++i);
	if (opt_G)
	    printf(count == 0 ? ",%d" : " %d",i+1);
	else
	    printf("%s ",Lat_CfName(&LI,i));
	++count;
    }
    else
    {
	PrintCompositionSeries(n->sub);
	PrintCompositionSeries(n->quot);
    }
}


/* ------------------------------------------------------------------
   WriteResult() - Write some information to stdout and create the
	.cfinfo file.
   ------------------------------------------------------------------ */

static void WriteResult(node_t *root)

{
    int i, k;

    MESSAGE(0,("\n\nChopping completed: %d different composition"
	" factors\n",LI.NCf));

    /* Write the cfinfo file
       --------------------- */
    MESSAGE(0,("Writing %s.cfinfo\n",LI.BaseName));
    for (i = 0; i < LI.NCf; ++i)
    {
    	LI.Cf[i].dim = irred[i]->dim;
	LI.Cf[i].num = irred[i]->num;
	LI.Cf[i].mult = irred[i]->mult;
	LI.Cf[i].spl = irred[i]->spl;
	LI.Cf[i].idword = irred[i]->idword;
	LI.Cf[i].idpol = irred[i]->idpol;
    }
    Lat_WriteInfo(&LI);


    /* Write composition factors
       ------------------------- */
    if (opt_G)
    {
	printf("MeatAxe.CompositionFactors := [\n");
    	for (i = 0; i < LI.NCf; ++i)
    	{
	    printf("  [ \"%s\", %ld, %ld ]",Lat_CfName(&LI,i),
	      irred[i]->mult,irred[i]->spl);
	    if (i < LI.NCf-1) printf(",");
	    printf("\n");
	}
	printf("\n];\n");
    }
    else if (MSG0)
    {
    	printf("\nName   Mult  SF  Fingerprint\n");
    	for (i = 0; i < LI.NCf; ++i)
    	{
	    printf("%-6s %4ld  %2ld  ",Lat_CfName(&LI,i),
	      irred[i]->mult,irred[i]->spl);
	    for (k = 0; k < MAXFP; ++k)
		printf("%d%s",irred[i]->fprint[k],k==MAXFP-1?"\n":",");
	}
    }

    /* Write the composition series
       ---------------------------- */
    if (opt_G)
    {
	printf("MeatAxe.CompositionSeries := [\n");
    	PrintCompositionSeries(root);
	printf("];\n");
    }
    else if (MSG0)
    {
	printf("\nAscending composition series:\n");
    	PrintCompositionSeries(root);
	printf("\n");
    }

    /* Write statistics
       ---------------- */
    if (MSG2)
    {
	printf("\nStatistics\n----------\n");
	printf("Saved vectors split: %4ld\n",stat_svsplit);
	printf("c(x) irreducible:    %4ld\n",stat_cpirred);
	printf("Normal split:        %4ld\n",stat_nssplit);
	printf("Dual split:          %4ld\n",stat_dlsplit);
	printf("Exceptional split:   %4ld\n",stat_exsplit);
	printf("Irreducible:         %4ld\n",stat_irred);
    }
}


/* --------------------------------------------------------------------------
   splitnode() - Split a node

   Arguments:
     <n>: Pointer to the node to split.
     <tr>: Indicates that it was a `dual split'

   Description:
     This function is called after a proper submodule has been found. It
     calculates the action of the generators on both submodule and quotient,
     and creates two new <note_t> structures for the two parts.
     The original module, <n>, is cleaned up and <Chop()> is called with
     both submodule and quotient.
   -------------------------------------------------------------------------- */

static void splitnode(node_t *n, int tr)

{	
    MatRep_t *sub = NULL, *quot = NULL;

    /* Split the module
       ---------------- */
    MESSAGE(0,("Split: Subspace=%d, Quotient=%d\n",
	n->subsp->Nor,n->dim - n->subsp->Nor));
    Split(n->subsp,tr ? n->TrRep : n->Rep,&sub,&quot);

    /* If it was a dual split, subspace and quotient have been
       calculated in the dual module. To get back to the original
       module, transpose again and exchange sub and quot.
       ---------------------------------------------------------- */
    if (tr)
    {	
	int i;
	Matrix_t *x, *y;
	for (i = 0; i < LI.NGen; ++i)
	{
	    x = MatTransposed(sub->Gen[i]);
	    MatFree(sub->Gen[i]);
	    y = MatTransposed(quot->Gen[i]);
	    MatFree(quot->Gen[i]);
	    sub->Gen[i] = y;
	    quot->Gen[i] = x;
	}
    }
    
    /* Make new nodes for subspace and quotient
       ---------------------------------------- */
    n->sub = CreateNode(sub,n);
    n->quot = CreateNode(quot,n);

    /* Project saved vectors on the quotient
       ------------------------------------- */
    if (!tr && n->nsp != NULL)
	n->quot->nsp = QProjection(n->subsp,n->nsp);

    /* Clean up
       -------- */
    CleanUpNode(n,1);

    /* Chop the subspace and quotient
       ------------------------------ */
    Chop(n->sub);
    Chop(n->quot);
}


/* ------------------------------------------------------------------
   extendbasis() - Find a vector in 'space' which is not in the
	(linear) span of 'basis'. 'basis' must be linearly
	independent.
   ------------------------------------------------------------------ */

static Matrix_t *extendbasis(Matrix_t *basis, Matrix_t *space)

{
    long i, j, piv;
    long dimb = basis->Nor;
    long dims = space->Nor;
    PTR tmp, x, y;
    FEL f;


    /* Concatenate basis and space
       --------------------------- */
    FfSetNoc(basis->Noc);
    tmp = FfAlloc(dimb + dims);
    memcpy(tmp,basis->Data,FfCurrentRowSize * dimb);
    /*x = FfGetPtr(tmp,dimb,FfNoc);*/
    x = (PTR)((char*)tmp + dimb * FfCurrentRowSize);
    memcpy(x,space->Data,FfCurrentRowSize * dims);

    /* Clean with basis
       ---------------- */
    for (i = 0, x = tmp; i < dimb; FfStepPtr(&x), ++i)
    {
	piv = FfFindPivot(x,&f);
	if (piv == -1) 
	    MTX_ERROR("extendbasis(): zero vector in basis");
	y = x;
	for (j = i+1; j < dimb+dims; ++j)
	{
	    FfStepPtr(&y);
	    FfAddMulRow(y,x,FfSub(FF_ZERO,FfDiv(FfExtract(y,piv),f)));
	}
    }


    /* Find the first non-zero row
       --------------------------- */
/*    x = FfGetPtr(tmp,dimb,FfNoc);*/
    x = (PTR)((char*)tmp + dimb * FfCurrentRowSize);
    for (j = 0; FfFindPivot(x,&f) == -1; ++j, FfStepPtr(&x));
    free(tmp);
    if (j > dims)
    {
	MTX_ERROR("extendbasis() failed");
	return NULL;
    }
    return MatCutRows(space,j,1);
}



/* ------------------------------------------------------------------
   checkspl() - Checks if a given representation's splitting field
   has degree [E:F] = dim(V), where V is a given subspace (usually
   the kernel of an algebra element).
   
   Returns 1 if [E:F]=dim(V) or zero otherwise.
   ------------------------------------------------------------------ */

#define MAXENDO 10	/* Max. dimension of endomorphism ring */


static int checkspl(const MatRep_t *rep, Matrix_t *nsp)

{
    Matrix_t *sb1, *sb2;	/* Standard bases */
    Matrix_t *g1[MAXGEN];	/* Generators in standard basis sb1 */
    Matrix_t *g2[MAXGEN];	/* Generators in standard basis sb2 */
    MatRep_t *endo;
    int i, result = 0;

    MESSAGE(2,("checkspl(): nsp=%d\n",nsp->Nor));
    /* Take the first vector from nsp and change to standard basis.
       ------------------------------------------------------------ */
    sb1 = SpinUp(nsp,rep,SF_FIRST|SF_CYCLIC|SF_STD,NULL,NULL);
    MTX_VERIFY(sb1 != NULL && sb1->Nor == sb1->Noc);
    ChangeBasisOLD(sb1,LI.NGen,(const Matrix_t **)rep->Gen,g1);
    endo = MrAlloc(0,NULL,0);

    sb2 = NULL;	/* Mark as unused */
    while (1)
    {  
	Matrix_t *v2, *subsp;

	/* Spin up v1 under all endomorphisms found so far. If this
	   yields the whole null-space, we know that the endomorphism
	   ring has at least dimension dim(nsp).
	   ---------------------------------------------------------- */
        MESSAGE(3,("1st spin-up: nendo=%d\n",endo->NGen));
	subsp = SpinUp(nsp,endo,SF_FIRST|SF_SUB,NULL,NULL);
	MTX_VERIFY(subsp != NULL);
	if (subsp->Nor == nsp->Nor)
	{
	    MatFree(subsp);
	    result = 1;	/* Successfull! */
	    break;
	}

	/* Take a vector which is not in span(v1)
	   and make again the standard basis.
	   --------------------------------------- */
        MESSAGE(3,("1st spin-up not successful\n"));
	if (sb2 != NULL)	/* Clean up first */
	{
	    int j;
	    MatFree(sb2);
	    for (j = 0; j < LI.NGen; ++j) 
		MatFree(g2[j]);
	}
	v2 = extendbasis(subsp,nsp);	/* Get vector */
	MatFree(subsp);
	sb2 = SpinUp(v2,rep,SF_FIRST|SF_CYCLIC|SF_STD,NULL,NULL);
	MTX_VERIFY(sb2 != NULL && sb2->Nor == sb2->Noc);
	MatFree(v2);
	ChangeBasisOLD(sb2,rep->NGen,(const Matrix_t **)rep->Gen,g2);

	/* Compare the two representations. If they are different,
	   we know that the splitting field degree must be smaller
	   than dim(nsp).
	   -------------------------------------------------------- */
	result = 1;
	for (i = 0; result != 0 && i < LI.NGen; ++i)
	{
	    if (MatCompare(g1[i],g2[i]) != 0)
	        result = 0;
	}
	if (result == 0) break;	/* Not successfull */

	/* They are identical, i.e., we have found an endomorphism.
	   Put it into the list and try the next vector.
	   -------------------------------------------------------- */
	MrAddGenerator(endo,MatMul(MatInverse(sb2),sb1),0);
    }

    /* Clean up
       -------- */
    MatFree(sb1); 
    for (i = 0; i < LI.NGen; ++i) MatFree(g1[i]);
    if (sb2 != NULL)
    {
	MatFree(sb2);
        for (i = 0; i < LI.NGen; ++i) MatFree(g2[i]);
    }
    MrFree(endo);
    MESSAGE(3,("checkspl(): result = %d\n",result));
    return result;
}



/* --------------------------------------------------------------------------
   FindIdWord() - Find an identifying word for a module
   
   Description:
     This function finds an identifying word for a given irreducible module. 
     This is a word in the generators with minimal nullity, i.e., the nullity 
     equals the splitting field degree [E:F]. The word is stored in n->idword, 
     the polynomial in n->idpol, and its null-space in n->nsp.

   Arguments:
     <n>: The module

   Return:
    0 = Ok, -1 = Error.
   -------------------------------------------------------------------------- */

static int FindIdWord(node_t *n)
{
    long i;
    int k;
    long count = 0;
    static FPoly_t *cpol = NULL;

    /* Main loop: Try all words
       ------------------------ */
    MESSAGE(1,("Searching idword, dim=%d\n",n->dim));
    for (i = 1; count <= MAXTRIES; ++i)
    {
	if (IsBadWord(i,n)) 
	    continue;

	/* Make the word and its characteristic polynomial
	   ----------------------------------------------- */
	MESSAGE(2,("  Word %ld  gcd=%ld\n",i,n->ggt));
	MakeWord(n,i);
	if (cpol != NULL) FpFree(cpol);
	cpol = CharPol(n->word);
	if (MSG3) FpPrint("  c(x)",cpol);
	for (k = 0; k < (int) cpol->NFactors; ++k)
	    n->ggt = gcd(n->ggt,cpol->Mult[k] * cpol->Factor[k]->Degree);

	/* Try all factors of c(x) with degree<=gcd
	   ---------------------------------------- */
	for (k = 0; k < (int) cpol->NFactors; ++k)
	{
	    if (cpol->Factor[k]->Degree > n->ggt) 
		continue;
	    ++count;
	    if (MSG3) PolPrint("  factor",cpol->Factor[k]);
	    InsertWord(n,cpol->Factor[k]);
	    MESSAGE(3,("  nsp=%d",n->nsp->Nor));
	    n->ggt = gcd(n->ggt,n->nsp->Nor);
	    if (n->nsp->Nor > n->ggt)
		continue;

	    if (checkspl(n->Rep,n->nsp))
	    {
		MESSAGE(2,("  idword = %ld\n",i));
		n->idword = i;
		SetInsert(goodwords,i);
		n->idpol = PolDup(cpol->Factor[k]);
		return 0;
	    }
	}
    }

    /* Failed...
       --------- */
    MTX_ERROR("FindIdWord() failed");
    MrSave(n->Rep,"CHOPFAILED");
    return -1;
}




/* ------------------------------------------------------------------
   newirred() - Check if a given irreducible module is already
	contained in the list of composition factors. If yes, return
	the index. If not, insert the new irreducible and return its
	index.
   ------------------------------------------------------------------ */

static void newirred(node_t *n)

{
    int i, k;
    Matrix_t *b;
    
    /* Check if the module is already in the list
       ------------------------------------------ */
    WgMakeFingerPrint(n->wg,n->fprint);
    for (i = 0; i < LI.NCf && n->dim >= irred[i]->dim; ++i)
    {
	/* Compare dimensions and fingerprints
	   ----------------------------------- */
	if (n->dim != irred[i]->dim ||
	    memcmp(n->fprint,irred[i]->fprint,sizeof(n->fprint)))
	    continue;

	if (IsIsomorphic(irred[i]->Rep,LI.Cf + i,n->Rep,NULL,0))
	{
	    ++irred[i]->mult;
	    n->num = irred[i]->num;
	    MESSAGE(0,("Irreducible (%s)\n",Lat_CfName(&LI,i)));
	    CleanUpNode(n,0);
	    return;
	}
    }

    /* It's a new irreducible!
       ----------------------- */
    if (LI.NCf >= LAT_MAXCF)
	MTX_ERROR("TOO MANY CONSTITUENTS");
    for (k = LI.NCf-1; k >= i; --k)
    {
	irred[k+1] = irred[k];
	LI.Cf[k+1] = LI.Cf[k];
    }
    irred[i] = n;
    ++LI.NCf;

    /* Set the fields
       -------------- */
    if (i == 0 || irred[i]->dim != irred[i-1]->dim)
	irred[i]->num = 0;
    else
	irred[i]->num = irred[i-1]->num + 1;
    irred[i]->mult = 1;
    LI.Cf[i].dim = irred[i]->dim;  /* cfname() needs this */
    LI.Cf[i].num = irred[i]->num;
    MESSAGE(0,("Irreducible (%s)\n",Lat_CfName(&LI,i)));

    /* Make idword and change to std basis 
       ----------------------------------- */
    MATFREE(n->nsp);
    if (n->idword == -1)
	FindIdWord(n);
    else
    {
	Matrix_t *m;
    	n->word = WgMakeWord(n->wg,n->idword);
    	m = MatInsert(n->word,irred[i]->idpol);
    	n->nsp = MatNullSpace__(m);
    }
    LI.Cf[i].idword = n->idword;
    LI.Cf[i].idpol = n->idpol;
    LI.Cf[i].spl = n->spl = n->nsp->Nor;
    b = SpinUp(n->nsp,n->Rep,SF_FIRST|SF_CYCLIC|SF_STD,NULL,NULL);
    MTX_VERIFY(b != NULL && b->Nor == b->Noc);
    ChangeBasisOLD(b,LI.NGen,(const Matrix_t **)n->Rep->Gen,n->Rep->Gen);
    MatFree(b);

    /* Write out the generators
       ------------------------ */
    if (irred[i]->spl > 1)
	MESSAGE(1,("Splitting field has degree %ld\n",irred[i]->spl));
    for (k = 0; k < LI.NGen; ++k)
    {
	char fn[200];
	sprintf(fn,"%s%s.%d",LI.BaseName,Lat_CfName(&LI,i),k+1);
   	MatSave(irred[i]->Rep->Gen[k],fn);
    }
    CleanUpNode(n,0);
}



/* ------------------------------------------------------------------
   SplitWithSavedVectors() - Try to split a module using saved
   vectors.
   Returns 1 on success, 0 otherwise.
   ------------------------------------------------------------------ */

static int SplitWithSavedVectors(node_t *n)

{
    Matrix_t *span;

    if (n->nsp == NULL || n->nsp->Nor == 0) 
	return 0;
    MESSAGE(1,("Trying saved vectors..."));

    if (n->subsp != NULL) 
	MatFree(n->subsp);
    n->subsp = NULL;
    span = SpinUp(n->nsp,n->Rep,SF_EACH|SF_SUB,NULL,NULL);
    if (span->Nor > 0 && span->Nor < span->Noc)
    {
	n->subsp = span;
	MESSAGE(1,("success\n"));
	++stat_svsplit;
	splitnode(n,0);
	return 1;
    }

    MatFree(span);
    MESSAGE(1,("failed\n"));
    return 0;
}


/* --------------------------------------------------------------------------
   polymap() - Calculate the image of a single vector -- or a set
   of vectors -- under p(A).
   -------------------------------------------------------------------------- */

Matrix_t *polymap(Matrix_t *v, Matrix_t *m, Poly_t *p)

{
    Matrix_t *result, *tmp;
    int i;

    result = MatAlloc(v->Field,v->Nor,v->Noc);
    tmp = MatDup(v);
    for (i = 0; i <= p->Degree; ++i)
    {
	FEL f = p->Data[i];
	PTR x = result->Data, y = tmp->Data;
	long i;
	for (i = v->Nor; i > 0; --i)
	{   
	    FfAddMulRow(x,y,f);
	    FfStepPtr(&x);
	    FfStepPtr(&y);
	}
	MatMul(tmp,m);
    }
    MatFree(tmp);
    return result;
}




/* --------------------------------------------------------------------------
   make_kern() - Calculate a vector in the null-space of p(A). p(x)
   is assumed to be a factor of f1(x), i.e., p(x) must occur in the
   first cyclic subspace.
   -------------------------------------------------------------------------- */

static void make_kern(node_t *n, Poly_t *p)

{
    Poly_t *cof, *f;
    Matrix_t *seed;

    f = PolDup(n->f1);
    cof = PolDivMod(f,p);
    MTX_VERIFY(f->Degree == -1);
    PolFree(f);

    if (n->nsp != NULL) 
	MatFree(n->nsp);
    seed = MatAlloc(FfOrder,1,n->dim);
    FfInsert(seed->Data,CharPolSeed,FF_ONE);
    n->nsp = polymap(seed,n->word,cof);
    MatFree(seed);
    PolFree(cof);
}


/* ------------------------------------------------------------------
   make_trkern() - Calculate a vector in the null-space of p(A^T). We
   consider only the first cyclic subspace, so make_trkern() may not
   find a vector.

   Returns 0 if a vector has been found, -1 else.
   ------------------------------------------------------------------ */

int make_trkern(node_t *n, Poly_t *p)

{
    Matrix_t *mt;
    Poly_t *pt, *cof;
    int result = -1;

    mt = MatTransposed(n->word);
    pt = CharPolFactor(mt);
    cof = PolDivMod(pt,p);
    if (pt->Degree == -1)		/* Glueck gehabt! */
    {
	Matrix_t *seed = MatAlloc(FfOrder,1,n->dim);
    	FfInsert(seed->Data,CharPolSeed,FF_ONE);
    	if (n->nsptr != NULL) 
	    MatFree(n->nsptr);
    	n->nsptr = polymap(seed,mt,cof);
	MatFree(seed);
	result = 0;
    }
    PolFree(cof);
    PolFree(pt);
    MatFree(mt);
    return result;
}


/* ------------------------------------------------------------------
   make_f1() - Calculate the characteristic polynomial on the first
   cyclic subspace.
   ------------------------------------------------------------------ */

static FPoly_t *make_f1(node_t *n)

{
    FPoly_t *cpol;
    int i, k;

    if (n->f2 != NULL) PolFree(n->f2);
    n->f2 = NULL;
    if (n->f1 != NULL) PolFree(n->f1);
    n->f1 = CharPolFactor(n->word);
    cpol = Factorization(n->f1);

    /* Sort the factors by ascending multiplicity
       ------------------------------------------ */
    for (i = 0; i < cpol->NFactors; ++i)
    {
    	for (k = i+1; k < (int)cpol->NFactors; ++k)
    	{
	    long val1 = cpol->Mult[i];
	    long val2 = cpol->Mult[k];
	    if (val1 > val2)
	    {
	    	Poly_t *p = cpol->Factor[i];
	    	long e = cpol->Mult[i];
	    	cpol->Factor[i] = cpol->Factor[k];
	    	cpol->Mult[i] = cpol->Mult[k];
	    	cpol->Factor[k] = p;
	    	cpol->Mult[k] = e;
	    }
    	}
    }
    if (MSG3) FpPrint("  f1(x)",cpol);
    return cpol;
}


/* --------------------------------------------------------------------------
   make_f2() - Complete the characteristic polynomial

   Description:
     This function completes the calculation of the characteristic
     polynomial. We assume that the first factor (i.e. the characteristic
     polynomial on the first cyclic subspace) has already be calculated.
     The remaining factors are stored in <f2>.
   --------------------------------------------------------------------------- */

static void make_f2(node_t *n)	/* Make f2 */
{
    Poly_t *f;

    if (n->f2 != NULL) 
	return;
    n->f2 = PolAlloc(n->f1->Field,0);
    while ((f = CharPolFactor(NULL)) != NULL)
    {
	PolMul(n->f2,f);
	PolFree(f);
    }
    if (MSG3)
    { 
	FPoly_t *x = Factorization(n->f2);
	FpPrint("  f2(x)",x);
	FpFree(x);
    }
}


/* --------------------------------------------------------------------------
   SplitWithNsp() - Try to split the module with the kernel in <n->nsp>.

   Return: 1 = success, 0 = failed
   -------------------------------------------------------------------------- */


static int SplitWithNsp(node_t *n)

{
    Matrix_t *sub;

    MATFREE(n->subsp);
    MESSAGE(3,("Trying to split with null-space..."));
    sub = SpinUp(n->nsp,n->Rep,SF_FIRST|SF_SUB,NULL,NULL);
    if (sub->Nor > 0 && sub->Nor < sub->Noc)
    {
	/* Module was split.
	   ----------------- */
	MESSAGE(3,("Success\n"));
	n->subsp = sub;
	++stat_nssplit;
	SetInsert(goodwords,n->wnum);
	splitnode(n,0);
	return 1;
    }
    MatFree(sub);
    MESSAGE(3,("Failed\n"));
    return 0;
}





static int PolMultiplicity(const Poly_t *factor, const Poly_t *pol)

{
    int mult = 0;
    Poly_t *f2 = PolDup(pol);
    int done = 0;
    while (!done)
    {
	Poly_t *q = PolDivMod(f2,factor);
	if (f2->Degree == -1)
	{
	    ++mult;
	    PolFree(f2);
	    f2 = q;
	}
	else
	{
	    PolFree(q);
	    done = 1;
	}
    }
    return mult;
}




/* ------------------------------------------------------------------
   try_poly() - Tries to split with one irreducible factor of c(x).
   
    Return: 1 = success, 0 = failed
   ------------------------------------------------------------------ */

static int try_poly(node_t *n, Poly_t *pol, long vfh)

{
    int chkirred = 1;
    Matrix_t *sub;
    int mult;

    /* Try to split.
       ------------- */
    make_kern(n,pol);
    if (SplitWithNsp(n))
	return 1;

    /* Now find out if we can prove irreducibility with Norton's
       criterion. We need that the factor we just tried occurs
       with multiplicity 1 in the characteristic polynomial, c(x).
       ----------------------------------------------------------- */
    if (n->f2 == NULL) 
	make_f2(n);
    mult = PolMultiplicity(pol,n->f2) + vfh;	/* Multiplicity in c(x) */

    /* If the multiplicity is not one, we can still prove irreducibility,
       but we must spin up every vector in the null-space of p(A). We do
       this only if the null-space is 'small' (as defined by -n).
       ------------------------------------------------------------------ */
    if (mult > 1)
    {
	if (mult * pol->Degree > opt_nullimit)
	{
	    chkirred = 0;
	    MESSAGE(3,("Cannot check for irreducibility, null-space=%d\n",
		mult * pol->Degree));
	}
	else
	{
	    /* Spin up all vectors */
	    MatFree(n->nsp);
	    n->nsp = MatNullSpace__(MatInsert(n->word,pol));
	    MATFREE(n->subsp);
	    MESSAGE(3,("2nd spin-up, null-space = %d\n",n->nsp->Nor));
	    sub = SpinUp(n->nsp,n->Rep,SF_MAKE|SF_SUB,NULL,NULL);
	    if (sub->Nor > 0 && sub->Nor < sub->Noc)
	    {
		/* Module was split */
		n->subsp = sub;
		++stat_nssplit;
		SetInsert(goodwords,n->wnum);
		splitnode(n,0);
		return 1;
	    }
	}
    }

    /* Not split, try dual split
       ------------------------- */
    if (make_trkern(n,pol) != 0)
    {
	MESSAGE(3,("No seed vector found, dual split skipped\n"));
	return 0;
    }
    MESSAGE(3,("Try dual split..."));
    if (n->TrRep == NULL)	/* Transpose generators */
	n->TrRep = MrTransposed(n->Rep);
    sub = SpinUp(n->nsptr,n->TrRep,SF_FIRST|SF_SUB,NULL,NULL);
    if (sub->Nor > 0 && sub->Nor < sub->Noc)
    {
	MESSAGE(3,("Success!\n"));
	n->subsp = sub;
	++stat_dlsplit;
	splitnode(n,1);
	return 1;
    }
    MatFree(sub);
    MESSAGE(3,("Failed\n"));
    if (!chkirred)
	return 0;

    /* The module is irreducible!
       -------------------------- */
    newirred(n);
    SetInsert(goodwords,n->wnum);
    ++stat_irred;
    return 1;
}



/* --------------------------------------------------------------------------
   try_ex_factor() - Try to split in exceptional cases

   Arguments:
     <n>: Pointer to the module
     <cp>: Characteristic polynomial
     <cpf>: Characteristic polynomial in factored form
     <factor>: Index of the factor to try

   Return:
     1 = success (module was split), 0 = failed

   Description:
     This function is called from <try_exceptional()> for each irreducible
     factor of the characteristic polynomial of degree >= 2.
   -------------------------------------------------------------------------- */

static int try_ex_factor(node_t *n, Poly_t *cp, FPoly_t *cpf, int factor)

{
    int k;
    long rndword;
    Poly_t *p, *q, *tmp, *gcd[3];
    Matrix_t *B, *iA, *v, *v1;

    if (MSG3) 
    {
	printf("Trying factor (");
	PolPrint(NULL,cpf->Factor[factor]);
	printf(")^%d\n",cpf->Mult[factor]);
    }

    /* Calculate p(x) = maximal power of the irred factor in c(x)
       ---------------------------------------------------------- */
    p = PolDup(cpf->Factor[factor]);
    for (k = 1; k < cpf->Mult[factor]; ++k)
	PolMul(p,cpf->Factor[factor]);

    /* Calculate the complement q(x) with q(x)p(x)=c(x)
       ------------------------------------------------ */
    tmp = PolDup(cp);
    q = PolDivMod(tmp,p);
    MTX_VERIFY(tmp->Degree == -1);
    PolFree(tmp);

    /* Calculate i(x):=b(x)q(x) with a(x)p(x)+b(x)q(x) = 1.
       ---------------------------------------------------- */
    PolGcdEx(p,q,gcd);
    MTX_VERIFY(gcd[0]->Degree == 0);
    PolMul(q,gcd[2]);

    /* Insert the word into i(x) and clean up polynomials.
       --------------------------------------------------- */
    iA = MatInsert(n->word,q);
    PolFree(gcd[0]);
    PolFree(gcd[1]);
    PolFree(gcd[2]);
    PolFree(p);
    PolFree(q);

    /* Choose a second random word, B, and calculate [A,i(A)Bi(A)]
       ----------------------------------------------------------- */
    rndword = n->wnum + MtxRandomInt(42);
    MESSAGE(3,("Choosing random word %ld\n",rndword));
    B = WgMakeWord(n->wg,rndword);
/* OLD
    MatMul(B,iA);
    MatMul(iA,B);
    MatFree(B);
    B = MatDup(n->word);
    MatMul(B,iA);
    MatMul(iA,n->word);
    MatMulScalar(iA,FfSub(FF_ZERO,FF_ONE));
    MatAdd(B,iA);
    MatFree(iA);
*/

    /* Select a random vector in the image of the commutator
       ----------------------------------------------------- */
    MESSAGE(3,("Choosing random vector\n"));
    v = MatAlloc(B->Field,1,B->Noc);
    for (k = 0; k < v->Noc; ++k)
	FfInsert(v->Data,k,FfFromInt(MtxRandomInt(FfOrder)));


    v1 = MatDup(v);
    MatMul(MatMul(MatMul(MatMul(v,n->word),iA),B),iA);
    MatMul(MatMul(MatMul(MatMul(v1,iA),B),iA),B);
    MatAddMul(v,v1,FfNeg(FF_ONE));
    MatFree(v1);
    MatFree(iA);
/* OLD
    MatMul(v,B);
*/
    MatFree(B);

    /* Try to split with this vector
       ----------------------------- */
    MESSAGE(3,("Spinning up:"));
    if (n->subsp != NULL) 
	MatFree(n->subsp);
    n->subsp = SpinUp(v,n->Rep,SF_FIRST|SF_SUB,NULL,NULL);
    MatFree(v);
    if (n->subsp->Nor > 0 && n->subsp->Nor < n->subsp->Noc)
    {
	MESSAGE(3,("SUCCESS!\n"));
	++stat_exsplit;
    	splitnode(n,0);
	return 1;
    }
    MESSAGE(3,("failed\n"));
    return 0;
}


/* --------------------------------------------------------------------------
   try_exceptional() - Try to split in exceptional cases

   Arguments:
     <n>: Pointer to the module

   Return:
     1 = success (module was split), 0 = failed

   Description:
     This function tries to split a module. It uses an algorithm developed by
     G. Ivanyos and K. Lux which is specifically designed for the `exceptional'
     cases where the standard methods fail.
   -------------------------------------------------------------------------- */

static int try_exceptional(node_t *n)

{
    Poly_t *cp;
    FPoly_t *cpf;
    int factor;
    int success = 0;
    MTX_VERIFY(n->f1 != NULL && n->f2 != NULL);

    MESSAGE(2,("Trying exceptional cases\n"));

    /* Calculate the complete characteristic polynomial c(x) 
       and its irreducible factors.
       ----------------------------------------------------- */
    cp = PolDup(n->f1);
    PolMul(cp,n->f2);
    cpf = Factorization(cp);

    /* Try all factors of degree >= 2
       ------------------------------ */
    for (factor = 0; factor < (int)cpf->NFactors && !success; ++factor)
    {
	if (cpf->Factor[factor]->Degree < 2)
	    continue;
	success = try_ex_factor(n,cp,cpf,factor);
    }

    FpFree(cpf);
    PolFree(cp);
    return success;
}




/* --------------------------------------------------------------------------
   ChopWithWord() - Chop a module using a given word. Returns 1 on
   success, 0 otherwise.
   -------------------------------------------------------------------------- */

static int ChopWithWord(node_t *n, long wn, int try_ex)

{
    long dlimit = opt_deglimit;		/* Limit on degree */
    int pi, done;
    FPoly_t *cpol;

    MESSAGE(2,("Next word is %ld (=%s)\n",wn,WgSymbolicName(n->wg,wn)));
    CharPolSeed = (CharPolSeed + 2) % n->dim;	/* BUG: '+2' for 2.3 compat. */
    MakeWord(n,wn);
    cpol = make_f1(n);		/* Make first part of c(x) */

    /* If c(x) is irreducible, then the module is irreducible
       ------------------------------------------------------ */
    if (cpol->Factor[0]->Degree == n->dim)
    {
	MESSAGE(2,("c(x) is irreducible\n"));
	newirred(n);
	++stat_cpirred;
	SetInsert(goodwords,wn);
	FpFree(cpol);
	return 1;
    }

    /* Try all factors of c(x)
       ----------------------- */
    for (done = pi = 0; !done && pi < (int) cpol->NFactors; ++pi)
    {    
	if (MSG3) {
	    printf("Next factor: ");
	    PolPrint(NULL,cpol->Factor[pi]);
	    printf("^%d\n",cpol->Mult[pi]);
	}
	if (dlimit > 0 && cpol->Factor[pi]->Degree > dlimit)
	{
	   MESSAGE(3,("deg > %ld -- discarded\n",dlimit));
	   continue;
	}
	done = try_poly(n,cpol->Factor[pi],cpol->Mult[pi]);
	MESSAGE(3,("try_poly()=%d\n",done));

    }
    FpFree(cpol);
    if (!done && try_ex)
        done = try_exceptional(n);
    return done;
}









/* --------------------------------------------------------------------------
   IsLinear() - Handle 1-dimensional modules. If the module is 1-dimensional,
   insert it into the list of irrecurcibles.

   Return: 0 on success, 1 if dim > 1.
   -------------------------------------------------------------------------- */

static int IsLinear(node_t *n)

{
    if (n->dim > 1)
	return 0;
    MESSAGE(1,("Dimension is one -- irreducible!\n"));
    newirred(n);
    ++stat_irred;
    return 1;
}



/* ------------------------------------------------------------------
   Chop() - Chop a module
   ------------------------------------------------------------------ */

static int Chop(node_t *n)

{
    long w;
    int count = 0;

    if (n == NULL)
    {
	MTX_ERROR1("node=NULL: %E",MTX_ERR_BADARG);
	return -1;
    }
    MESSAGE(0,("Chop: Dim=%d\n",n->dim));

    /* Some special cases first: 1-dimensional, saved vectors.
       ------------------------------------------------------- */
    if (IsLinear(n))
	return 0;
    if (SplitWithSavedVectors(n))
	return 0;

    /* Try words. First, we try all known 'good' words (w=0,-1,-2,...),
       then all words (w=1,2,...). Known 'bad' words are always skipped.
       ----------------------------------------------------------------- */
    MESSAGE(1,("Trying words\n"));
    for (count = 0, w = 0; count < 10000000; ++count)
    {
	long i;
	if ((w <= 0) && (int)-w < goodwords->Size)
	    i = goodwords->Data[(int) -(w--)];
	else 
	{
	    if (w <= 0)
		w = i = firstword;	/* Start with word 1 */
	    else
		i = ++w;
	    if (SetContains(goodwords,w))  /* Do not try again */
		continue;
	}
	if (IsBadWord(i,n)) 
	    continue;
	if (ChopWithWord(n,i,count > 10))
	    break;
	SetInsert(n->badwords,i);
    }
    return 0;
}



static int Init(int argc, const char **argv)

{
    goodwords = SetAlloc();
    memset(&LI,0,sizeof(Lat_Info));
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;

    /* Parse command line.
       ------------------- */
    opt_G = AppGetOption(App,"-G --gap");
    opt_i = AppGetOption(App,"-i --read-cfinfo");
    firstword = AppGetIntOption(App,"-s",1,1,100000);
    LI.NGen = AppGetIntOption(App,"-g",2,1,MAXGEN);
    opt_deglimit = AppGetIntOption(App,"-d",-1,-1,100);
    opt_nullimit = AppGetIntOption(App,"-n",FfOrder > 10 ? 2 : 3,1,20);
    if (AppGetArguments(App,1,1) != 1)
	return -1;
    strcpy(LI.BaseName,App->ArgV[0]);
    return 0;
}


static void Cleanup()

{
    if (App != NULL)
	AppFree(App);
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, const char **argv)

{
    
    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }

    /* Chop the module
       --------------- */
    MESSAGE(0,("*** CHOP MODULE ***\n\n"));
    if (CreateRoot() != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }
    Chop(root);

    /* Write result, clean up, and exit
       -------------------------------- */
    WriteResult(root);

    Cleanup();
    return 0;
}

/**
@page prog_chop chop - Find Irreducible Constituents

@section chop_syntax Command Line
<pre>
chop [@em Options] [-Gi] [-g @em NGen] [-s @em Word] [-n @em MaxNul] [-d @em MaxDeg] @em Name
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par @em -G --gap
Produce GAP output. Implies -Q.

@par @em -i --read-cfinfo
Read @em Name.cfinfo, if it exists, to determine the number of generators.

@par -g @em NGen
Set the number of generators. Default is two generators, but see -i.

@par -s @em Word
Start with the given word number instead of 1.

@par -n @em MaxNul
Set limit on nullity. Only null-spaces with a dimension less than or equal to
@em MaxNul are searched completely.

@par -d @em MaxDeg
Set limit on degrees of polynomials.

@par @em Name
Name of the module to chop.

@section chop_inp Input Files
@par @em Name.1, @em Name.2, ...
Action of the generators on the module.

@section chop_out Output Files

@par @em Name.cfinfo
Constituent information file

@par @em Name.X.1, @em Name.X.2, ...
Action of the generators on the constituent X.

@section chop_desc Description
This program calculates the irreducible constituents of a given matrix representation.
The representing matrices of the generators are read from input files, see "input Files"
above. Unless a different number of generators has been specified with -g, two
generators are expected. However, if the "-i" option is used, and the file @em Name.cfinfo
exists, @b chop takes the number of generators from this file and ignores the "-g" option.

For each composition factor @b chop writes the action of the generators to @em CFName.1,
@em CFName.2, ... @em CFName is the name of the composition factor, which is constructed
by appending the dimension and a letter to the module name. For example, "X10a.1"
is the action of the first generator on the first composition factor of dimension
10 of the module X. If a second, inequivalent composition factor of dimension 10 
was found, it would be named `X10b' and so on.
@b chop also creates the file @em Name.cfinfo' containing a list of all composition factors.
This file is used by subsequent programs such as @ref prog_pwkond "pwkond".

@section chop_impl Implementation Details
@b chop repeatedly splits a module into submodule and quotient until it arrives at the
irreducible constituents. Thus, it finds a composition series. The program assumes that the
algebra generated by the input matrices contains the unit matrix.

In order to split a given module or to prove its irreducibility the algorithm needs an element
of the algebra with a non-trivial but low-dimensional kernel. Such elements are searched by
taking linear combinations of certain products of the generators ("words").
See the description of the @ref prog_zmw "zmw" program for more details on the word generator.
By default, @b chop tries all words in the order defined by the word generator. The "-s" option
may be used to make @b chop start with a word different from 1.

For each word A generated in this way, the program calculates its characteristic polynomial and
examines the irreducible factors.  If p(x) is an irreducible factor, p(A) has a non-trivial
kernel.  Then, one vector of the kernel is chosen and the smallest submodule containing this
vector is calculated.
If the vector spans a proper submodule, the action of the generators on this submodule as well
as on the quotient are calculated and the same procedure is applied recursively to both
submodule and quotient.

To avoid expensive matrix multiplications in the calculation of p(A),
there is a limit on the degree of p(x). This limit can be set with
the "-d" option and defaults to 5.

If a module cannot be split by the program, it may be irreducible.
In order to prove this, @b chop uses Norton's criterion. This requires,
however, to find an algebra element with a small kernel, because up
to scalar multiples each vector in the kernel must be examined to see
whether it spins up to the whole module. For this reason a "nullity
threshold" m is maintained by the program. Initially, m is set to 3 or
to the value given in the "-n" option.
Each algebra element that has a nullity less then or
equal to m is used for the Norton test.

In some cases the algorithm described is not able to split the module
although it is reducible. These exceptional cases are treated with
an alternative strategy described in @ref LI98 "[LI98]".

Algebra elements with trivial kernel are useless for the algorithm, so
an attempt is made to avoid unnecessary computation of such elements.
Once an element is known to have a trivial kernel on a given module
M, the program will mark it as invertible and ignore it for all
constituents of M.

If a constituent is irreducible but not absolutely irreducible, the
nullity of any element in the algebra will be a multiple of [E:F],
where F is the ground field and E the splitting field.
This situation is recognized by calculating the greatest common divisor d
of all nullities which occur during the search.
In order to prove that the splitting field degree is equal to
@em d, the following method is used: Take a word with nullity @em d
and two vectors v1, v2 in its null-space. Use these vectors as
seeds for a standard basis algorithm. If the resulting representations
are different, [E:F] is less than @em d, and the word is discarded.
Otherwise, the linear map
which transforms one standard basis into the other is an endomorphism
@em e of the module. If v1, under the action of e, spins up to
the whole null space, then [E:F]=@em d. Otherwise, take a third vector
not in the span and repeat the procedure above. Again, this yields an
endomorphism, or it turns out that [E:F]<@em d.
These steps are repeated until a word with nullity [E:F] is found.
*/
