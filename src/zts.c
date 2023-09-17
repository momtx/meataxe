////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tensor spin
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */




static MtxApplicationInfo_t AppInfo = { 
"zts", "Tensor Spin", 
"\n"
"SYNTAX\n"
"    zts [<Options>] <M> <N> <Seed> [<Sub>]\n"
"\n"
"ARGUMENTS\n"
"    <M> ..................... Left representation\n"
"    <N> ..................... Right representation\n"
"    <Seed> .................. Seed vector(s)\n"
"    <Sub> ................... Invariant subspace\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -g <#Gens> .............. Set number of generators (default: 2)\n"
"    -n, --no-action ......... Output only <Sub>, do not calculate <Sub>.1, ...\n"
"\n"
"FILES\n"
"    <M>.{1,2...} ............ I Generators on left representation\n"
"    <N>.{1,2...} ............ I Generators on right representation\n"
"    <Seed> .................. I Seed vectors\n"
"    <Sub>.................... O Invariant subspace\n"
"    <Sub>.{1,2...} .......... O Action on the invariant subspace\n"
};

static MtxApplication_t *App = NULL;


/**
 ** Internal data structure.
 **/

typedef struct
{
    int Row;
    int Col;
    FEL Mark;
} tPivotEntry;


static const char *MName = "?";	    /* Name of left factor (M) */
static const char *NName = "?";	    /* Name of right factor (N) */
static const char *SeedName = "?";  /* Name of seed vector file */
static const char *SubName = "?";   /* Name of result file */
static int NGen = 2;		    /* Number of generators */
static Matrix_t **GenM = NULL;	    /* Generators on M (transposed) */
static Matrix_t **GenN = NULL;	    /* Generators on N */
static Matrix_t **Basis = NULL;	    /* Basis of invariant subspace */
static tPivotEntry *Piv = NULL;	    /* Has alwas the same size as Basis */
static int Dim = 0;		    /* Number of M-Vectors in <Basis> */
static int Src = 0;		    /* Index of vector to map next (0-based) */
static int MaxDim = 0;		    /* Capacity of <Basis> */ 
static Matrix_t *Seed;		    /* Seed vectors */
static int NoAction = 0;	    /* -n: Only subspace, no actions */
static int TpDim = 0;		    /* Dimension of the tensor product */

static int ReadFiles()
{
    int i;
    char fn[200];
    Matrix_t *m;

    GenM = NALLOC(Matrix_t *,NGen);
    GenN = NALLOC(Matrix_t *,NGen);
    for (i = 0; i < NGen; ++i)
    {
	sprintf(fn,"%s.%d",MName,i+1);
	m = matLoad(fn);
	GenM[i] = matTransposed(m);
	matFree(m);
	sprintf(fn,"%s.%d",NName,i+1);
	if ((GenN[i] = matLoad(fn)) == NULL)
	    return -1;
    }
    Seed = matLoad(SeedName);
    if (Seed  == NULL)
	return -1;
    TpDim = GenM[0]->Nor * GenN[0]->Nor;
    MESSAGE(1,("Tensor product has dimension %d*%d=%d\n",
	GenM[0]->Nor,GenN[0]->Nor,TpDim));
    return 0;
}


static void VecToMat(PTR vec, Matrix_t *m)
{
    int row, pos;
    PTR rowptr;

    for (pos=row=0, rowptr=m->Data; row < m->Nor; ++row,ffStepPtr(&rowptr,m->Noc))
    {
	int col;
	ffMulRow(rowptr,FF_ZERO, m->Noc);
	for (col = 0; col < m->Noc; ++col, ++pos)
	{
	    FEL f = ffExtract(vec,pos);
	    if (f != FF_ZERO)
		ffInsert(rowptr,col,f);
	}
    }
}

static void matToVec(PTR vec, const Matrix_t *m)
{
    int row, pos;
    PTR rowptr;

    ffMulRow(vec,FF_ZERO, TpDim);
    for (pos=row=0, rowptr=m->Data; row < m->Nor; ++row,ffStepPtr(&rowptr,m->Noc))
    {
	int col;
	for (col = 0; col < m->Noc; ++col, ++pos)
	{
	    FEL f = ffExtract(rowptr,col);
	    if (f != FF_ZERO)
		ffInsert(vec,pos,f);
	}
    }
}


static int FindPivot(Matrix_t *m, tPivotEntry *piv)
{
    int row;
    PTR rowptr;

    for (row = 0, rowptr=m->Data; row < m->Nor; ++row, ffStepPtr(&rowptr,m->Noc))
    {
	int col;
	FEL f;
	col = ffFindPivot(rowptr,&f,m->Noc);
	if (col != MTX_NVAL)
	{
	    piv->Row = row;
	    piv->Col = col;
	    piv->Mark = f;
	    return 0;
	}
    }
    return -1;			/* Not found, zero matrix */
}


static void Clean(Matrix_t *mat, const Matrix_t **basis, 
    const tPivotEntry *piv, int dim)
{
    int i;

    MTX_ASSERT(dim == 0 || (mat->Noc == basis[0]->Noc && mat->Nor == basis[0]->Nor));
    for (i = 0; i < dim; ++i)
    {
	PTR x;
	FEL f;
	x = matGetPtr(mat,piv[i].Row);
	f = ffNeg(ffDiv(ffExtract(x,piv[i].Col),piv[i].Mark));
	matAddMul(mat,basis[i],f);
    }
}



static void Clean2(Matrix_t *mat, const Matrix_t **basis, 
    const tPivotEntry *piv, int dim, PTR op)

{
    int i;

    MTX_ASSERT(mat->Noc == basis[0]->Noc && mat->Nor == basis[0]->Nor);
    ffMulRow(op,FF_ZERO, dim);
    for (i = 0; i < dim; ++i)
    {
	PTR x;
	FEL f;
	x = matGetPtr(mat,piv[i].Row);
	f = ffDiv(ffExtract(x,piv[i].Col),piv[i].Mark);
	matAddMul(mat,basis[i],ffNeg(f));
	ffInsert(op,i,f);
    }
}


static Matrix_t *Map(Matrix_t *src, int gen)
{
    Matrix_t *a;

    MTX_ASSERT(gen >= 0 && gen < NGen);
    a = matDup(GenM[gen]);
    matMul(a,src);
    matMul(a,GenN[gen]);
    return a;
}



static void CleanAndAppend(Matrix_t *mat)
{
    tPivotEntry newpiv;

    Clean(mat,(const Matrix_t **)Basis,Piv,Dim);
    if (FindPivot(mat,&newpiv) == 0)
    {
	if (Dim >= MaxDim)
	{
	    MaxDim += 50;
	    Basis = NREALLOC(Basis,Matrix_t *,MaxDim);
	    Piv = NREALLOC(Piv,tPivotEntry,MaxDim);
	}
	Piv[Dim] = newpiv;
	Basis[Dim] = mat;
	++Dim;
	if (Dim % 100 == 0)
	    MESSAGE(2,("Dimension=%d (%d%%)\n",Dim,Src*100/Dim));
    }
    else
	matFree(mat);
}



static void SpinUpMatrix(Matrix_t *seed)
{
    int gen = 0;
    Matrix_t *newvec = seed;

    Src = Dim;
    while (1)
    {
        CleanAndAppend(newvec);
	if (Src >= Dim)
	    break;		/* Done! */
	newvec = Map(Basis[Src],gen);
	if (++gen >= NGen)
	{
	    gen = 0;
	    ++Src;
	}
    }
}



static void Spinup()
{
    int i;
    PTR vec;

    Basis = NULL;
    Piv = NULL;
    MaxDim = 0;
    Dim = 0;
    vec = Seed->Data;
    for (i = 1; i <= Seed->Nor; ++i)
    {
	Matrix_t *seed;
	MESSAGE(1,("Spinning up seed vector %d\n",i));
	seed = matAlloc(ffOrder,GenM[0]->Nor,GenN[0]->Nor);
	VecToMat(vec,seed);
	SpinUpMatrix(seed);		    /* <Spinup()> eats <seed>! */
	ffStepPtr(&vec,Seed->Noc);
	if (i < Seed->Nor)
	    MESSAGE(1,("Dimension = %d\n",Dim));
    }
    MESSAGE(0,("Subspace has dimension %d\n",Dim));
}


static void WriteSubspace()
{
    int i;
    MtxFile_t *f;
    PTR row;

    MESSAGE(1,("Writing subspace to %s\n",SubName));
    row = ffAlloc(1, TpDim);
    f = mfCreate(SubName,Seed->Field,Dim,TpDim);
    for (i = 0; i < Dim; ++i)
    {
	matToVec(row,Basis[i]);
	mfWriteRows(f,row,1,TpDim);
    }
    mfClose(f);
    ffFree(row);
}


static void CalculateAction1(int gen, const char *file_name)
{
    MtxFile_t *f;
    PTR rowptr;
    int i;

    MESSAGE(1,("Writing generator to %s\n",file_name));
    f = mfCreate(file_name,Seed->Field,Dim,Dim);
    rowptr = ffAlloc(1, Dim);
    for (i = 0; i < Dim; ++i)
    {
	Matrix_t *image = Map(Basis[i],gen);
	Clean2(image,(const Matrix_t **)Basis,Piv,Dim,rowptr);
	matFree(image);
	mfWriteRows(f,rowptr,1,Dim);
    }
    ffFree(rowptr);
    mfClose(f);
}




static void CalculateAction()
{
    int i;
    MESSAGE(0,("Calculating action of generators on subspace\n"));
    for (i = 0; i < NGen; ++i)
    {
	char fn[200];
	sprintf(fn,"%s.%d",SubName,i+1);
	CalculateAction1(i,fn);
    }
}



static int Init(int argc, char **argv)
{
    /* Process command line options.
       ----------------------------- */
    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    NGen = appGetIntOption(App,"-g",2,1,1000);
    NoAction = appGetOption(App,"-n --no-action");

    /* Procerss command line arguments.
       -------------------------------- */
    if (appGetArguments(App,3,4) < 0)
	return -1;
    MName = App->ArgV[0];
    NName = App->ArgV[1];
    SeedName = App->ArgV[2];
    SubName = App->ArgC >= 4 ? App->ArgV[3] : NULL;
    return 0;
}



static void Cleanup()
{
    appFree(App);
}




int main(int argc, char **argv)
{
    if (Init(argc, argv) != 0)
	return 1;
    if (ReadFiles() != 0)
    {
	mtxAbort(MTX_HERE,"Cannot read input file(s)");
	return 1;
    }
    Spinup();
    if (SubName != NULL)
    {
	WriteSubspace();
	if (!NoAction)
	    CalculateAction();
    }
    Cleanup();
    return 0;
}



/**
@page prog_zts zts - Tensor Split


@section zts_syntax Command Line
<pre>
zts @em Options [-g @em NumGen] @em M @em N @em Seed [@em Sub]
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par -g
    Set number of generators (default is 2).

@par -n --no-action
    Do not calculate the action of the generators on the invariant
    subspace. Output only the subspace.

@par @em M
    First representation (left factor).

@par @em N
    Second representation (right factor).

@par @em Seed
    Seed vector(s).

@par @em Sub
    Invariant subspace. Also used as basename for the action
    on the invariant subspace.


@section zts_inp Input Files

@par \em M.1, \em M.2, ...
Generator action on left module. Unless changed with -g,
two generators are read.

@par @em N.1, @em N.2, ...
Generator action on right module

@par @em Seed
Seed vector(s).

@section zts_out Output Files

@par @em Sub
Basis of the invariant subspace.

@par @em Sub.1, @em Sub.2, ...
Generator action on the invariant subspace.

@see
- @ref prog_tuc
- @ref prog_zsp

@section zts_desc Description
This program is similar to @ref prog_zsp "zsp", but it works on the tensor
product of two modules, M⊗N. @b zts spins up one or more vectors, and optionally
calculates a matrix representation corresponding to the invariant
subspace. The program does not use the matrix representation of the
generators on M⊗N, which would be too large in many cases.
This program is used, for example, to spin up vectors that have
been uncondensed with @ref prog_tuc "tuc".

The action of the generators on both @em M and @em N must be given as square matrices,
see "Input Files" above.
You can use the -g option to specify the number of generators.
The default is two generators.

Seed vectors are read from @em Seed.
They must be given with respect to the lexicographically ordered
basis explained below.
    
If the @em Sub argument is given, ZTS writes a basis of the
invariant subspace to @em Sub, calculates the action of the
generators on the invariant subspace, and writes it to
@em Sub.1, @em Sub.2,...

@section zts_impl Implementation Details
Let \f$B=(b_1,\ldots,b_m)\f$ be a basis of \e M,
\f$C=(c_1,\ldots,c_n)\f$ a basis of \e N,
and denote by \f$B\otimes C\f$ the lexicographically ordered basis
\f$(b_1\otimes c_1, b_1\otimes c_2, \dots b_m\otimes c_n)\f$.
For \f$v\in M\otimes N\f$,
the coordinate row \f$m(v,B\otimes C)\f$ has \f$mn\f$ entries
which can be arranged as a \f$m\times n\f$ matrix (top to bottom, left to
right).  Let \f$M(B,v,C)\f$ denote this matrix. Then
\f[
	M(B,va,C) = m(B,a|_M,B)^{\rm tr}M(B,v,C)m(C,a|_N,C)
	\quad\mbox{for all }
	a\in A,\ v\in M\otimes N
\f]
Using this relation,
we can calculate the image of any vector \f$v\in M\otimes N\f$ under an algebra element \f$a\f$,
and thus spin up a vector without using the matrix representation of \f$a\f$
on \f$v\in M\otimes N\f$.

*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
