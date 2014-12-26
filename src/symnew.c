#include "meataxe.h"
#include <stdlib.h>
#include <string.h>

MTX_DEFINE_FILE_INFO

#define MAX_DEGREE 5
#define MAX_NPERMS 10

/// @todo **/

typedef struct
{
    int Factor;
    int Perm[MAX_DEGREE];
} FPerm_t;

/// @todo **/

typedef struct
{
    int Degree;
    int Nominator;
    int NPerms;
    FPerm_t Fp[MAX_NPERMS];
} Symmetrizer_t;

Symmetrizer_t E3 = {
    3, /* Degree */
    6, /* Nominator */
    6, { {  1, { 0,1,2} }, { 1, { 1,2,0} }, { 1, { 2,0,1} },
         { -1, { 1,0,2} }, {-1, { 0,2,1} }, {-1, { 2,1,0} } }
};

Symmetrizer_t M3 = {
    3, /* Degree */
    3, /* Nominator */
    3, { {  2, { 0,1,2} }, { -1, { 1,2,0} }, { -1, { 2,0,1} } }
};

Symmetrizer_t S2 = {
    2, /* Degree */
    2, /* Nominator */
    2, { {  1, { 0,1} }, { 1, { 1,0} } } 
};

int Degree;
int VDim, WDim;


/*-----------------------------------------------------------*/
int OmegaLen;
FEL OmegaFactor[MAX_NPERMS];
const int *OmegaPerm[MAX_NPERMS];


static void MakeOmega(const Symmetrizer_t *s)
{
    FEL f;
    int i;
    OmegaLen = s->NPerms;
    f = FfFromInt(s->Nominator);
    for (i = 0; i < OmegaLen; ++i)
    {
	OmegaFactor[i] = FfDiv(FfFromInt(s->Fp[i].Factor),f);
	OmegaPerm[i] = s->Fp[i].Perm;
    }
}

/*-----------------------------------------------------------*/

static void NumToTuple(int num, int *tuple)
{
    int i;
    MTX_ASSERT(num >= 0 && num < WDim);
    for (i = Degree - 1; i >= 0; --i)
    {
	tuple[i] = num % VDim;
	num /= VDim;
    }
}

static int TupleToNum(const int *base)
{
    int i;
    int num = 0;
    for (i = Degree - 1; i >= 0; --i)
    {
	num = num * VDim + *base;
	++base;
    }
    MTX_ASSERT(num >= 0 && num < WDim);
    return num;
}
/*-----------------------------------------------------------*/

/// @todo **/
typedef struct entry { int Num; FEL Coeff; } SvEntry_t;
/// @todo **/
typedef struct {
    int Size;
    int MaxSize;
    SvEntry_t *V;
} SvVector_t;

int SDim = 0;
SvVector_t **SBasis = NULL;


SvVector_t *SvAlloc(int size)
{
    SvVector_t *x = ALLOC(SvVector_t);
    memset(x,0,sizeof(*x));
    x->Size = 0;
    x->MaxSize = size;
    x->V = NALLOC(SvEntry_t,size);
    memset(x->V,0,size * sizeof(SvEntry_t));
    return x;
}

int SvFree(SvVector_t *vec)
{
    SysFree(vec->V);
    SysFree(vec);
    return 0;
}

static int SvFind(const SvVector_t *vec, int n, int *pos)
{
    int i;
    int diff = 1;
    for (i = 0; i < vec->Size && (diff = n - vec->V[i].Num) > 0; ++i)
	;
    *pos = i;
    return diff;
}


void SvAddEntry(SvVector_t *vec, int n, FEL f)
{
    int pos;
    if (f == FF_ZERO)
	return;
    if (SvFind(vec,n,&pos) == 0)
    {
	vec->V[pos].Coeff = FfAdd(vec->V[pos].Coeff,f);
	if (vec->V[pos].Coeff == FF_ZERO)
	{
	    int k;
	    for (k = pos + 1; k < vec->Size; ++k)
		vec->V[k-1] = vec->V[k];
	    --vec->Size;
	}
    }
    else
    {
	int k;
	if (vec->Size >= vec->MaxSize)
	{
	    vec->MaxSize = vec->Size + 1;
	    vec->V = NREALLOC(vec->V,SvEntry_t,vec->MaxSize);
	}
	for (k = vec->Size - 1; k >= pos; --k)
	    vec->V[k+1] = vec->V[k];
	vec->V[pos].Num = n;
	vec->V[pos].Coeff = f;
	++vec->Size;
    }
} 

void SvPrint(const SvVector_t *vec)
{
    int i;

    for (i = 0; i < vec->Size; ++i)
    {
	int tuple[MAX_NPERMS];
	int k;
	NumToTuple(vec->V[i].Num,tuple);
	if (i > 0)
	    printf("+");
	printf("%dv[",FfToInt(vec->V[i].Coeff));
	for (k = 0; k < Degree; ++k)
	    printf(k == 0 ? "%d" : ",%d",tuple[k]);
	printf("]");
    }
}

void SvAddMul(SvVector_t *vec, const SvVector_t *b, FEL f)
{
    int i;
    if (f == FF_ZERO)
	return;
    for (i = 0; i < b->Size; ++i)
	SvAddEntry(vec,b->V[i].Num,FfMul(b->V[i].Coeff,f));
}


void SvClean2(SvVector_t *vec, const SvVector_t **basis, int dim, PTR op)
{
    int i;
    if (op != NULL)
	FfMulRow(op,FF_ZERO);
    for (i = 0; i < dim; ++i)
    {
	FEL f;
	int num = basis[i]->V[0].Num;	/* Pivot 'column' */
	int pos;
    	if (SvFind(vec,num,&pos) != 0)
	    continue;
	f =  FfDiv(vec->V[pos].Coeff,basis[i]->V[0].Coeff);
	SvAddMul(vec,basis[i],FfNeg(f));
    	if (op != NULL)
	    FfInsert(op,i,f);
    }	
}



void SvClean(SvVector_t *vec, const SvVector_t **basis, int dim)
{
    SvClean2(vec,basis,dim,NULL);
}


/*-----------------------------------------------------------*/

static void MapBasisVector(int n, SvVector_t *x)
{
    int tuple[MAX_DEGREE];
    int imgtuple[MAX_DEGREE];
    int i;

    NumToTuple(n,tuple);
    x->Size = 0;    

    for (i = 0; i < OmegaLen; ++i)
    {
	int k, img_n;
	for (k = 0; k < Degree; ++k)
	    imgtuple[OmegaPerm[i][k]] = tuple[k];
	img_n = TupleToNum(imgtuple);
	SvAddEntry(x,img_n,OmegaFactor[i]);
    }
}

#if 0
static void PrintRow(PTR row)
{
    int i;
    const char *fmt = (FfOrder < 10) ? "%d" : " %d";
    for (i = 0; i < FfNoc; ++i)
	printf(fmt,FfExtract(row,i));
}	
#endif


static void MakeSBasis()
{
    int num;
    SvVector_t *v = SvAlloc(OmegaLen);
    SBasis = NALLOC(SvVector_t *,WDim);   /* Worst case */
    
    SDim = 0;
    for (num = 0; num < WDim; ++num)
    {
	MapBasisVector(num,v);
	SvClean(v,(const SvVector_t **)SBasis,SDim);
	if (v->Size > 0)
	{
	    /* Add new basis vector */
	    if (MSG2)
	    {
	    	printf("SBasis[%d]=",SDim);
	    	SvPrint(v);
	    	printf("\n");
	    }
	    SBasis[SDim++] = v;
	    v = SvAlloc(OmegaLen);
	}
    }
    SvFree(v);
    /* OPT: free unused memory in <SBasis> */
    MESSAGE(0,("Dim(S) = %d\n",SDim));
}


static void BuildTables(const Symmetrizer_t *s, int dim_v)
{
    int i;

    Degree = s->Degree;
    VDim = dim_v;
    for (i = Degree - 1, WDim = VDim; i > 0; --i)
	WDim *= VDim;
    MESSAGE(0,("Dim(V) = %d\n",VDim));
    MESSAGE(0,("Degree = %d\n",Degree));
    MESSAGE(0,("Dim(W) = %d\n",WDim));

    MakeOmega(s);
    MakeSBasis();
}



MtxFile_t *OutputFile;

static void MapTuple(int t[], int s[], int start, FEL f, PTR mat,
	SvVector_t *result)
{
    int i;
    PTR row = FfGetPtr(mat,t[start]);
#if 0
    printf("Start=%d Row[%d]=",start,t[start]); PrintRow(row); printf("\n");
#endif
    for (i = 0; i < VDim; ++i)
    {
	FEL g = FfExtract(row,i);
	if (g == FF_ZERO)
	    continue;
	s[start] = i;
    	if (start < Degree - 1)
	    MapTuple(t,s,start+1,FfMul(f,g),mat,result);
	else
	{
	    int n = TupleToNum(s);
	    SvAddEntry(result,n,FfMul(f,g));
	}
    }
}


static void MapVector(int n, PTR mat, SvVector_t *result)
{
    int i;
    result->Size = 0;
    FfSetNoc(VDim);

    for (i = 0; i < SBasis[n]->Size; ++i)
    {
	int tuple[MAX_NPERMS];
	int imgtuple[MAX_NPERMS];
    	NumToTuple(SBasis[n]->V[i].Num,tuple);
	MapTuple(tuple,imgtuple,0,SBasis[n]->V[i].Coeff,mat,result);
    }
}

static void CalculateSAction(PTR mat)
{
    int n;
    PTR img;
    SvVector_t *v = SvAlloc(100);
    img = FfAlloc(1);
    for (n = 0; n < SDim; ++n)
    {
	MapVector(n,mat,v);
    	FfSetNoc(SDim);
	fprintf(stderr,".");
	SvClean2(v,(const SvVector_t **)SBasis,SDim,img);
	MfWriteRows(OutputFile,img,1);
    }
    SysFree(img);
}







static MtxApplicationInfo_t AppInfo = { 
"zsy", "Symmetrized Tensor Product",
"SYNTAX\n"
"    zsy " MTX_COMMON_OPTIONS_SYNTAX " [-G] <Mode> <Inp> <Out>\n"
"\n"
"ARGUMENTS\n"
"    <Mode> .................. Symmetrization mode: e2, e3, e4, s2, or m3\n"
"    <Inp> ................... Input matrix\n"
"    <Out> ................... Output matrix\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... GAP output (implies -Q)\n"
};

static MtxApplication_t *App = NULL;


static int opt_G = 0;		/* GAP output */
static const char *iname, *oname;	/* File names */
static enum {M_E2,M_E3,M_E4,M_S2,M_M3} mode;
static long fl;			/* Field */
static long nor, noc;		/* Input sizes */
static PTR m1;			/* Input */
static PTR *row;		/* Pointers to the rows of the input matrix */






/* ------------------------------------------------------------------
   init() - Initialize everything, read input files, allocate memory
   ------------------------------------------------------------------ */

static int Prepare()

{
    long i;
    MtxFile_t *f;

    /* Open the input file and check the header.
       ----------------------------------------- */
    if ((f = MfOpen(iname)) == NULL)
    {
	MTX_ERROR("Error opening input file");
	return 1;
    }
    if (f->Field < 2)
    {
	MTX_ERROR2("%s: %E",iname,MTX_ERR_NOTMATRIX);
	return 1;
    }
    if (f->Nor != f->Noc)
    {
	MTX_ERROR2("%s: %E",iname,MTX_ERR_NOTSQUARE);
	return 1;
    }
    fl = f->Field;
    nor = f->Nor;
    noc = f->Noc;

	FfSetField(fl); 
	FfSetNoc(noc);
	m1 = FfAlloc(nor);
	MfReadRows(f,m1,nor);

	/* Set up pointers to the rows of the input matrix
   	   ----------------------------------------------- */
	row = NALLOC(PTR,nor);
	if (row == NULL)
	    return -1;
	for (i = 0; i < nor; ++i)
	{	
	    row[i] = m1;
	    FfStepPtr(&m1);
	}
    MfClose(f);

    /*BuildTables(&M3,nor);
    BuildTables(&S2,nor);
      */
    BuildTables(&E3,nor);


    /* Allocate the output buffer
       -------------------------- */
    MESSAGE(0,("Output is %d x %d\n",SDim,SDim));
    FfSetNoc(SDim);
    OutputFile = MfCreate(oname,fl,SDim,SDim);


    return 0;
}



static int Init(int argc, const char **argv)

{
    const char *arg3;

    /* Process command line options.
       ----------------------------- */
    App = AppAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;
    opt_G = AppGetOption(App,"-G --gap");
    if (opt_G)
	MtxMessageLevel = -100;

    /* Process arguments.
       ------------------ */
    if (AppGetArguments(App,3,3) < 0)
	return -1;
    iname = App->ArgV[1];
    oname = App->ArgV[2];
    arg3 = App->ArgV[0];
    if (!strcmp(arg3,"e2")) mode = M_E2;
    else if (!strcmp(arg3,"e3")) mode = M_E3;
    else if (!strcmp(arg3,"e4")) mode = M_E4;
    else if (!strcmp(arg3,"s2")) mode = M_S2;
    else if (!strcmp(arg3,"m3")) mode = M_M3;
    else 
    {
	MTX_ERROR1("Unknown mode '%s'",arg3);
	return -1;
    }
    return 0;
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
    if (Prepare() != 0)
	return 1;
    CalculateSAction(row[0]);
    MfClose(OutputFile);
    AppFree(App);
    return EXIT_OK;
}


