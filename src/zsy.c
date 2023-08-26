////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Symmetric tensor product.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <stdlib.h>
#include <string.h>




/* ------------------------------------------------------------------
   Global Data
   ------------------------------------------------------------------ */


static MtxApplicationInfo_t AppInfo = { 
"zsy", "Symmetrized Tensor Product",
"SYNTAX\n"
"    zsy " MTX_COMMON_OPTIONS_SYNTAX " [-G] <Mode> <Inp> <Out>\n"
"\n"
"ARGUMENTS\n"
"    <Mode> .................. Symmetrization mode: e2, e3, e4, or s2\n"
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
static FILE *ofile = NULL;
static enum {M_E2,M_E3,M_E4,M_S2} mode;
static long fl;			/* Field */
static long nor, noc;		/* Input sizes */
static long nor2, noc2;		/* Output sizes */
static PTR m1;			/* Input */
static PTR m2;			/* One row of the output matrix */
static PTR *row;		/* Pointers to the rows of the input matrix */
static Perm_t *PermIn, *PermOut;






/* ------------------------------------------------------------------
   init() - Initialize everything, read input files, allocate memory
   ------------------------------------------------------------------ */

static int Prepare()

{
    long i;
    MtxFile_t *f;

    /* Open the input file and check the header.
       ----------------------------------------- */
    if ((f = mfOpen(iname)) == NULL)
    {
	mtxAbort(MTX_HERE,"Error opening input file");
	return 1;
    }
    if (mode != M_E2 && mode != M_S2 && mode != M_E3 && f->Field < 2)
    {
	mtxAbort(MTX_HERE,"%s: %s",iname,MTX_ERR_NOTMATRIX);
	return 1;
    }
    fl = f->Field;
    nor = f->Nor;
    noc = f->Noc;

    if (fl >= 2)	/* Matrix */
    {
	ffSetField(fl); 
	ffSetNoc(noc);
	m1 = ffAlloc(nor, noc);
	mfReadRows(f,m1,nor);

	/* Set up pointers to the rows of the input matrix
   	   ----------------------------------------------- */
	row = NALLOC(PTR,nor);
	if (row == NULL)
	    return -1;
	for (i = 0; i < nor; ++i)
	{	
	    row[i] = m1;
	    ffStepPtr(&m1, noc);
	}
    }
    else		/* Permutation */
    {	
	sysFseek(f->File,0);
	PermIn = permRead(f->File);
	if (PermIn == NULL) 
	    return -1;
    }
    mfClose(f);

    /* Calculate nor2 and noc2
       ----------------------- */
    switch (mode)
    {	case M_S2:
	    nor2 = (nor * (nor+1)) / 2;
	    noc2 = (noc * (noc+1)) / 2;
	    break;
	case M_E2:
	    nor2 = (nor * (nor-1)) / 2;
	    noc2 = (noc * (noc-1)) / 2;
	    break;
	case M_E3:
	    nor2 = nor*(nor-1)/2*(nor-2)/3;
	    noc2 = noc*(noc-1)/2*(noc-2)/3;
	    break;
	case M_E4:
	    nor2 = nor*(nor-1)/2*(nor-2)/3*(nor-3)/4;
	    noc2 = noc*(noc-1)/2*(noc-2)/3*(noc-3)/4;
	    break;
	default:
	    mtxAbort(MTX_HERE,"Unknown mode %d",(int)mode);
	    return -1;
    }
    if (fl > 0 && (nor2 < 0 || noc2 <= 0))
    {
	mtxAbort(MTX_HERE,"%s: Matrix too small",iname);
	return -1;
    }


    // Prepare output buffer and file.
    if (fl >= 2)
    {
	MESSAGE(0,("Output is %ld x %ld\n",nor2,noc2));
	ffSetNoc(noc2);
	if ((m2 = ffAlloc(1, noc2)) == NULL)
	    return -1;
	ofile = ffWriteHeader(oname,fl,nor2,(fl >= 2) ? noc2 : 1);
	if (ofile == NULL)
	{
	    mtxAbort(MTX_HERE,"Error creating output file");
	    return -1;
	}
    }
    else
    {	
	MESSAGE(0,("Output has degree %ld\n",nor2));
	fflush(stdout);
	PermOut = permAlloc(nor2);
	if (PermOut == NULL) 
	    return -1;
    }


    return 0;
}



/* ------------------------------------------------------------------
   zs2() - Symmetric square
   ------------------------------------------------------------------ */

static void zs2()
{   long i1, i2, j1, j2, j3;
    FEL f11, f12, f21, f22;
    FEL w1,w2,f1,f2;

    MESSAGE(1,("Mode S2, part 1\n"));
    for (i1 = 0; i1 < nor - 1; ++i1)
    {	
	for (i2 = i1 + 1; i2 < nor; ++i2)
	{   
	    ffMulRow(m2,FF_ZERO);
	    j3 = 0;
	    for (j1 = 0; j1 < noc - 1; ++ j1)
	    {	
		f11 = ffExtract(row[i1],j1);
		f21 = ffExtract(row[i2],j1);
		for (j2 = j1+1; j2 < noc; ++j2)
		{   
		    f12 = ffExtract(row[i1],j2);
		    f22 = ffExtract(row[i2],j2);
		    w1 = ffMul(f11,f22);
		    w2 = ffMul(f12,f21);
		    ffInsert(m2,j3,ffAdd(w1,w2));
		    ++j3;
		}
	    }
	    for (j2 = 0; j2 < noc; ++j2)
	    {	f1 = ffExtract(row[i1],j2);
		f2 = ffExtract(row[i2],j2);
		ffInsert(m2,j3,ffMul(f1,f2));
		++j3;
	    }
            MTX_ASSERT(ffNoc == noc2,);
	    ffWriteRows(ofile,m2,1, noc2);
	}
    }

    MESSAGE(1,("Mode S2, part 2\n"));
    for (i1 = 0; i1 < nor; ++i1)
    {	
	j3 = 0;
	ffMulRow(m2,FF_ZERO);
	for (j1 = 0; j1 < noc-1; ++j1)
	{   
	    f1 = ffExtract(row[i1],j1);
	    for (j2 = j1+1; j2 < noc; ++j2)
	    {	
		f2 = ffExtract(row[i1],j2);
		w2 = ffMul(f1,f2);
		ffInsert(m2,j3,ffAdd(w2,w2));
		++j3;
	    }
	}
	for (j2 = 0; j2 < noc; ++j2)
	{   
	    f1 = ffExtract(row[i1],j2);
	    ffInsert(m2,j3,ffMul(f1,f1));
	    ++j3;
	}
        MTX_ASSERT(ffNoc == noc2,);
	ffWriteRows(ofile,m2,1, noc2);
    }
}



/* ------------------------------------------------------------------
   zs2p() - Antisymmetric square (permutations)
   ------------------------------------------------------------------ */

static int maps2(int i, int k)

{	
    if (i <= k)
	return ((k*(k+1))/2 + i);
    else
	return ((i*(i+1))/2 + k);
}



static void zs2p()

{
    long *p1 = PermIn->Data;
    long *p2 = PermOut->Data;
    int i;

    for (i = 0; i < nor; ++i)
    {	
	int k;
	for (k = 0; k <= i; ++k)
	    p2[maps2(i,k)] = maps2(p1[i],p1[k]);
    }
    permSave(PermOut,oname);
}



/* ------------------------------------------------------------------
   ze2p() - Antisymmetric square (permutations)
   ------------------------------------------------------------------ */

static int mape2(int i, int k)

{
    if (i < k)
	return (k*(k-1))/2 + i;
    else
	return (i*(i-1))/2 + k;
}



static void ze2p()

{
    long *p1 = PermIn->Data;
    long *p2 = PermOut->Data;
    int i;

    for (i = 0; i < nor; ++i)
    {	
	int k;
	for (k = 0; k < i; ++k)
	    p2[mape2(i,k)] = mape2(p1[i],p1[k]);
    }
    permSave(PermOut,oname);
}




/* ------------------------------------------------------------------
   ze2() - Antisymmetric square
   ------------------------------------------------------------------ */

static void ze2()

{
    int i1, i2, j1, j2, j3;
    FEL f12, f22, w1, w2, w3;

    for (i1 = 0; i1 < nor-1; ++i1)
    {	
        MESSAGE(1,("i1 = %d\n",i1)); 
	for (i2 = i1+1; i2 < nor; ++i2)
	{
            MESSAGE(2,("i2 = %d\n",i2)); 
	    ffMulRow(m2,FF_ZERO);
	    j3 = 0;
	    for (j1 = 0; j1 < noc-1; ++j1)
	    {	register FEL f11, f21;
		f11 = ffExtract(row[i1],j1);
		f21 = ffExtract(row[i2],j1);
		for (j2 = j1+1; j2 < noc; ++j2)
		{
		    f12 = ffExtract(row[i1],j2);
		    f22 = ffExtract(row[i2],j2);
		    w1 = ffMul(f11,f22);
		    w2 = ffMul(f12,f21);
		    w3 = ffSub(w1,w2);
		    ffInsert(m2,j3,w3);
		    ++j3;
		}
	    }
            MTX_ASSERT(ffNoc == noc2,);
	    ffWriteRows(ofile,m2,1, noc2);
	}
    }
}


/* ------------------------------------------------------------------
   ze3() - Antisymmetric cube
   ------------------------------------------------------------------ */

static void ze3()

{   FEL f11,f12,f13,f21,f22,f23,f31,f32,f33;
    FEL e,g12,g13,g23;
    int i1, i2, i3, j1, j2, j3, jins;

    for (i1 = 0; i1 < nor-2; ++i1)
    {	
   	MESSAGE(1,("i1 = %d\n",i1)); 
	for (i2 = i1+1; i2 < nor-1; ++i2)
	{   
   	    MESSAGE(2,("i2 = %d\n",i2)); 
	    for (i3 = i2+1; i3 < nor; ++i3)
	    {	
   	       	MESSAGE(3,("i3 = %d\n",i3)); 
		ffMulRow(m2,FF_ZERO);
		jins = 0;
		for (j1 = 0; j1 < noc-2; ++j1)
		{   
		    f11 = ffExtract(row[i1],j1);
		    f21 = ffExtract(row[i2],j1);
		    f31 = ffExtract(row[i3],j1);
		    for (j2 = j1+1; j2 < noc-1; ++j2)
		    {	
			f12 = ffExtract(row[i1],j2);
			f22 = ffExtract(row[i2],j2);
			f32 = ffExtract(row[i3],j2);
			g12 = ffSub(ffMul(f11,f22),ffMul(f21,f12));
			g13 = ffSub(ffMul(f31,f12),ffMul(f11,f32));
			g23 = ffSub(ffMul(f21,f32),ffMul(f31,f22));
			for (j3 = j2+1; j3 < noc; ++j3)
			{   f13 =ffExtract(row[i1],j3);
			    f23 =ffExtract(row[i2],j3);
			    f33 =ffExtract(row[i3],j3);
			    e = ffAdd(ffAdd(ffMul(g12,f33),ffMul(g13,f23)),
			    	ffMul(g23,f13));
			    ffInsert(m2,jins,e);
			    ++jins;
			}
		    }
		}
                MTX_ASSERT(ffNoc == noc2,);
		ffWriteRows(ofile,m2,1l, noc2);
	    }
	}
    }
}


/* ------------------------------------------------------------------
   ze3p() - Antisymmetric cube (permutations)
   ------------------------------------------------------------------ */

#define SWAP(x,y) {tmp=x; x=y; y=tmp;}

static long mape3(long i, long k, long l)

{	
    register long tmp;
    if (i < k) SWAP(i,k);
    if (i < l) SWAP(i,l);
    if (k < l) SWAP(k,l);
    return (i*(i-1)/2*(i-2)/3 + k*(k-1)/2 + l);
}



static void ze3p()

{
    long *p1 = PermIn->Data;
    long *p2 = PermOut->Data;
    int i, k, l;
    for (i = 2; i < nor; ++i)
    {	
	for (k = 1; k < i; ++k)
	{
	    for (l = 0; l < k; ++l)
		p2[mape3(i,k,l)] = mape3(p1[i],p1[k],p1[l]);
	}
    }
    permSave(PermOut,oname);
}



/* ------------------------------------------------------------------
   ze4() - Antisymmetric fourth power
   ------------------------------------------------------------------ */

static void ze4()

{
   FEL f11,f12,f13,f14,f21,f22,f23,f24,f31,f32,f33,f34,f41,f42,f43,f44;
   FEL e,g12,g13,g14,g23,g24,g34,g123,g124,g134,g234;
   int i1, i2, i3, i4, j1, j2, j3, j4, jins;

   for (i1 = 0; i1 < nor-3; ++i1)
   {	
      MESSAGE(1,("i1 = %d\n",i1)); 
      for (i2 = i1+1; i2 < nor-2; ++i2)
      {   
         MESSAGE(2,("i2 = %d\n",i2)); 
         for (i3 = i2+1; i3 < nor-1; ++i3)
         {   
            MESSAGE(3,("i3 = %d\n",i3)); 
            for (i4 = i3+1; i4 < nor; ++i4)
            {   
               ffMulRow(m2,FF_ZERO);
               jins = 0;
               for (j1 = 0; j1 < noc-3; ++j1)
               {   
                  f11 = ffExtract(row[i1],j1);
                  f21 = ffExtract(row[i2],j1);
                  f31 = ffExtract(row[i3],j1);
                  f41 = ffExtract(row[i4],j1);

                  for (j2 = j1+1; j2 < noc-2; ++j2)
                  {   f12 = ffExtract(row[i1],j2);
                     f22 = ffExtract(row[i2],j2);
                     f32 = ffExtract(row[i3],j2);
                     f42 = ffExtract(row[i4],j2);

                     g12 = ffSub(ffMul(f11,f22),ffMul(f21,f12));
                     g13 = ffSub(ffMul(f11,f32),ffMul(f31,f12));
                     g14 = ffSub(ffMul(f11,f42),ffMul(f41,f12));
                     g23 = ffSub(ffMul(f21,f32),ffMul(f31,f22));
                     g24 = ffSub(ffMul(f21,f42),ffMul(f41,f22));
                     g34 = ffSub(ffMul(f31,f42),ffMul(f41,f32));

                     for (j3 = j2+1; j3 < noc-1; ++j3)
                     {   f13 =ffExtract(row[i1],j3);
                        f23 =ffExtract(row[i2],j3);
                        f33 =ffExtract(row[i3],j3);
                        f43 =ffExtract(row[i4],j3);

                        g123 = ffMul(f13,g23);
                        g123 = ffSub(g123,ffMul(f23,g13));
                        g123 = ffAdd(g123,ffMul(f33,g12));
                        g124 = ffMul(f13,g24);
                        g124 = ffSub(g124,ffMul(f23,g14));
                        g124 = ffAdd(g124,ffMul(f43,g12));
                        g134 = ffMul(f13,g34);
                        g134 = ffSub(g134,ffMul(f33,g14));
                        g134 = ffAdd(g134,ffMul(f43,g13));
                        g234 = ffMul(f23,g34);
                        g234 = ffSub(g234,ffMul(f33,g24));
                        g234 = ffAdd(g234,ffMul(f43,g23));

                        for (j4 = j3+1; j4 < noc; ++j4)
                        {   f14 =ffExtract(row[i1],j4);
                           f24 =ffExtract(row[i2],j4);
                           f34 =ffExtract(row[i3],j4);
                           f44 =ffExtract(row[i4],j4);

                           e = ffMul(f24,g134);
                           e = ffSub(e,ffMul(f14,g234));
                           e = ffAdd(e,ffMul(f44,g123));
                           e = ffSub(e,ffMul(f34,g124));
                           ffInsert(m2,jins,e);
                           ++jins;
                        }
                     }
                  }
               }
               MTX_ASSERT(ffNoc == noc2,);
               ffWriteRows(ofile,m2,1, noc2);
            }
         }
      }
   }
}





static int Init(int argc, char **argv)

{
    const char *arg3;

    /* Process command line options.
       ----------------------------- */
    App = appAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;
    opt_G = appGetOption(App,"-G --gap");
    if (opt_G)
	MtxMessageLevel = -100;

    /* Process arguments.
       ------------------ */
    if (appGetArguments(App,3,3) < 0)
	return -1;
    iname = App->ArgV[1];
    oname = App->ArgV[2];
    arg3 = App->ArgV[0];
    if (!strcmp(arg3,"e2")) mode = M_E2;
    else if (!strcmp(arg3,"e3")) mode = M_E3;
    else if (!strcmp(arg3,"e4")) mode = M_E4;
    else if (!strcmp(arg3,"s2")) mode = M_S2;
    else 
    {
	mtxAbort(MTX_HERE,"Unknown mode '%s'",arg3);
	return -1;
    }
    return 0;
}








/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char **argv)

{

    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }
    if (Prepare() != 0)
	return 1;
    switch (mode)
    {	
	case M_S2: fl >= 2 ? zs2() : zs2p(); break;
	case M_E2: fl >= 2 ? ze2() : ze2p(); break;
	case M_E3: fl >= 2 ? ze3() : ze3p(); break;
	case M_E4: ze4(); break;
	default:
	    mtxAbort(MTX_HERE,"Unknown mode %d",(int)mode);
	    return 1;
	    break;
    }
    if (ofile != NULL)
	fclose(ofile);
    appFree(App);
    return EXIT_OK;
}




/**
@page prog_zsy zsy - Symmetrized Tensor Product

@section zsy_syntax Command Line
<pre>
zsy [@em Options] [-G] @em Mode @em Inp @em Out
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -G
  GAP output.
@par @em Mode
  Symmetrization mode: "e2", "s2", "e3" or "e4".
@par @em Inp
  Input matrix.
@par @em Out
  Result matrix.

@section zsy_inp Input Files
@par @em Mat
  Input matrix or permutation.

@section zsy_out Output Files
@par @em Result
  Result matrix.

@section zsy_desc Description
This program reads a matrix or permutation, calculates its symmetrized tensor
product according to @em Mode, and writes out the result.

The @em Mode argument specifies the tensor product to be taken
and the kind of symmetrization to be performed. Currently there are
4 Modes available:
- "s2" is the symmetric tensor square. The output has size
  n(n+1)/2 (For matrices, number of lines, for permutations,
  degree).
- "e2" is the antisymmetric tensor square. The output has size
  n(n-1)/2.
- "e3" is the antisymmetric tensor cube. The output has size
  n(n-1)(n-2)/6.
- "e4" is the antisymmetric fourth power. The output has size
  n(n-1)(n-2)(n-3)/24.

Since the typical application of @b zsy is to generate new representations from 
existing ones, it will usually be used with square matrices. However,
the input is not required to be square.


@subsection zsy_perms Permutations
Currently, only modes s2, e2 and e3 are available for permutations.
The result gives the operation of the input permutation on unordered
pairs (e2, s2) or triples (e3) of points.
More precisely, if the given permutation operates on 1...n, then:
- s2 is the operation on (i,k) with 1≤i≤k≤n.
- e2 is the operation on (i,k) with 1≤i<k≤n.
- e3 is the operation on (i,k,l) with 1≤i<k<l≤n.

In the output, pairs and triples are numbered lexicographically.
For example, E2 uses the following order:
(1,2), (1,3), (2,3), (1,4), ...
Notice that the symmetric square is never transitive but 
decomposes into the diagonal and the antisymmetric square.
Here are some examples:
<pre>
p     = (1 5 4 3 2)
e2(p) = (1 7 10 6 3)(2 8 4 9 5)
s2(p) = (1 15 10 6 3)(2 11 14 9 5)(7 14 8 4 12)
e3(p) = (1 5 8 10 4)(2 6 9 3 7)
</pre>


@subsection mats Matrices
The r-th exterior power (modes e2, e3, e4) has as its entries the determinants of
r times r submatrices of the input. Rows and columns are ordered lexicographically,
which is equivalent to taking the following basis in the tensor product:
@par e2
	v<sub>i</sub> ∧ v<sub>j</sub> with 1≤i<j≤n
@par e3
	v<sub>i</sub> ∧ v<sub>j</sub> ∧ v<sub>k</sub> with 1≤i<j<k≤n
@par e4
	v<sub>i</sub> ∧ v<sub>j</sub> ∧ v<sub>k</sub>∧v<sub>l</sub> with 1≤i<j<k<l≤n

The basis vectors are ordered lexicographically, for example (e2):
v<sub>1</sub>∧v<sub>2</sub>, v<sub>1</sub>∧v<sub>3</sub>, ... v<sub>1</sub>∧v<sub>n</sub>,
v<sub>2</sub>∧v<sub>3</sub>, v<sub>2</sub>∧v<sub>4</sub>, ... v<sub>3</sub>∧v<sub>n</sub>,
... v<sub>n-1</sub>∧v<sub>n</sub>.

The symmetric square of a matrix with r rows and c columns is a
matrix with r(r+1)/2 rows and c(c+1)/2 columns, with entries
given by the formulae
@f[
\begin{array}{c|cc}
	  &  c(c-1)/2  & c  \\ \hline
r*(r-1)/2 & ad+bc      & ac \\
	r &  2ab       & a^2
\end{array}
@f]
where the upper left is the r(r-1)/2 by c(c-1)/2 matrix of permanents.
The program orders both the rows and the columns in lexicographical order, i.e.
v<sub>1</sub>v<sub>2</sub>, v<sub>1</sub>v<sub>3</sub>, ... v<sub>1</sub>v<sub>n</sub>,
v<sub>2</sub>v<sub>3</sub>, v<sub>2</sub>v<sub>4</sub>, ...
v<sub>2</sub>v<sub>n</sub>, v<sub>3</sub>v<sub>4</sub>, ... v_{n-1}v<sub>n</sub>,
v<sub>1</sub>v<sub>1</sub>, v<sub>2</sub>v<sub>2</sub>, ...
v<sub>n</sub>v<sub>n</sub>
with the assumption that v<sub>i</sub>v<sub>j</sub> = v<sub>j</sub>v<sub>i</sub>,
i.e. the action is on quadratic polynomials.

The symmetric square is, in general, irreducible except in characteristic 2.
In that case there is a copy of the Frobenius square
as an invariant submodule, as can be seen from the 2ab in the above
formulae. Invariant subspaces in characteristic 2 correspond to special
groups (i.e.\ groups of the form 2<sup>n</sup>×2<sup>m</sup>) on which the group
given acts on the quotient 2<sup>n</sup>.

Here are some examples:
<pre>
     (1 2 1 3)    (1 2 1 3 6 2)
  E2 (0 1 2 1) =  (0 1 0 2 0 4)     (mod 7)
     (1 2 2 3)    (6 5 6 5 1 4)

     (1 0 2 0 2)   (1 0 1 4 0 1 0 0 0 0)
  E3 (1 1 2 1 2) = (1 4 3 4 0 3 2 1 3 4)     (mod 5)
     (3 3 2 3 2)   (1 2 2 3 2 3 1 3 4 2)
     (1 2 3 1 0)   (4 0 4 0 2 0 4 2 1 3)

                   (1  2  1  5  5  7  0  2  2  3)
     (1 2 1 3)     (4  3  6  6 12  9  1  4  2  9)
  S2 (0 1 2 1)  =  (1  2  1  6  5  8  0  2  4  3)   (mod 13)
     (1 2 2 3)     (4  2  6  4 12  6  1  4  1  9)
                   (0  0  0  4  2  4  0  1  4  1)
                   (4  4  6  8 12 12  1  4  4  9)
</pre>




@section zsy_impl Implementation Details
If the input file contains more than one permutation, only the
first permutation is read in and processed.

If the input is a matrix, the whole input matrix and one row of the
result must fit into memory. In case of permutations both the input
and the result must fit into memory.

**/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
