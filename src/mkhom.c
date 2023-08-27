////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculate homomorphisms.
// Original program written by Magdolna Szőke. 
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

static MtxApplicationInfo_t AppInfo = { 
"mkhom", "Calculate homomorphisms",
"SYNTAX\n"
"    mkhom [-ts] [-r|-l] [-b] [-H <Dim>] <M> <N> <Hom>\n"
"\n"
"ARGUMENTS\n"
"    <M> ..................... First representation\n"
"    <N> ..................... Second representation\n"
"    <Hom> ................... Homomorhisms from <M> to <N>\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -t ...................... Calculate generators for <M> in spinning basis\n"
"    -s ...................... When <M>=<N>, give endomorphisms in spinning basis\n"
"    -r|-l ................... When <M>=<N>, find a generating set of End(M), and\n"
"                              calculate the left (-l) or right (-r) regular\n"
"                              representation.\n"
"    -b ...................... For big endorings, with -r, save memory.\n"
"    -H <Dim> ................ If the radical is given, <Dim> is the dimension of\n"
"\n"
"FILES\n"
"    <M>.{1,2...} ............ I  Generators in representation <M>.\n"
"    <N>.{1,2...} ............ I  Generators in representation <N>.\n"
"    <M>.cfinfo .............. I  Constituent info file for <M>.\n"
"    <N>.cfinfo .............. I  Constituent info file for <N>.\n"
"    <M>.rad ................. I  Generators for the head of <M> (with -H).\n"
"    <M><Cf>.k ............... I  Uncondense matrix, produced by PWKOND.\n"
"    <M>.std.................. O  The spinning basis for <M>.\n"
"    <Hom>.{1,2,...} ......... O  A k-basis of Hom(<M>,<N>).\n"
"    <M>.std.{1,2,...} ....... O  Generators in spinning basis (with -t).\n"
};

static MtxApplication_t *App = NULL;
static int standard = 0;                /* -t: Make standard generators */
static int hominstd = 0;                /* -s: Endomorphisms in std. basis */
static char reg = '?';                        /* 'l' or 'r' (with -r) */
static int big = 0;                        /* Argument of -b */
static int hd = 0;                        /* Head dimension (from -H) */
static const char *MName = NULL;        /* First argument */
static const char *NName = NULL;        /* Second argument */
static const char *HomName = NULL;        /* Third argument */
static int comp = 0;                        /* strcmp(MName,NName) */
static Lat_Info MInfo;                        /* Constituent info */
static MatRep_t *MRep = NULL;                /* First representation */
static MatRep_t *NRep = NULL;                /* Second representation */
static int dimM;                                /* Dimension of <M> */
static int partdim = 0;

static PTR basis = NULL;
static PTR space = NULL;
static int *piv = NULL;
static long *op = NULL;
static PTR *stdgen = NULL;
static long **stdtab = NULL;
static int *dims = NULL;
static Matrix_t *rad = NULL;
static Matrix_t ***posimages;        /* a basis of the possible images of the
                                                                   basis pieces */
static int *kerdim;
static Matrix_t *esys;
static int *esyspiv;
static int maxnumMgens = 0;
static int numMgens = 0;

/* TODO: ersetzen durch FfCleanRowAndRepeat */

static void myzcleanrow(PTR row, PTR matrix, PTR matrix2, int nor, int noc, int *piv)

{
    long i;
    PTR x, y, row2;

    row2 = ffGetPtr(matrix2,nor, noc);
    for (i = 0, x = matrix, y = matrix2; i < nor; 
         ++i, ffStepPtr(&x,noc), ffStepPtr(&y,noc))
    {
        FEL f = ffExtract(row, piv[i]);
        if (f == FF_ZERO) 
            continue;
        f = ffNeg(ffDiv(f, ffExtract(x, piv[i])));
        ffAddMulRow(row, x, f, noc);
        ffAddMulRow(row2, y, f, noc);
    }
}


/* ------------------------------------------------------------------
   zgensbasis() - Spin up canonically (spinning basis).

   <seed> ist ein Zeiger auf Seedvektoren, Der seedcount - 1-ste
   wird benutzt. <gen> sind die Erzeuger (gen[0]...gen[ngen-1]).

   <space> und <basis> muessen vom Aufrufer allokiert werden. Beide
   muessen gross genug fuer eine quadratrische Matrix sein. In <basis>
   wird die Spinningbasis abgelegt.

   <piv> ist die Pivot-Tabelle, die ebenfalls vom Aufrufer allokiert
   werden muss (mindestans dim + 1 Elemente!).

   In <op_table> wird die Definition der Spinningbasisvektoren abgelegt.
   (2 * (dim+1) long integers). Falls diese Information nicht
   benoetigt wird, kann op_table=NULL gesetzt werden.

   das Stueck, das der <seedcount>-ste Vektor ueber den schon existierenden
   <partdim> dimensionalen Modul erzeugt (der von den ersten seedcount - 1
   Vektoren aufgespannt ist) wird berechnet.

   <stdgen> sind die Standard erzeuger (stueckweise aufgebaut),
   in <stdtab> steht, welche Operationen zu ueberpruefen sind um zu entscheiden,
   ob ein gewisser Modul mit aehnlicher Basis zu dem gegebenen Modul isomorph
   ist.
   <stdgen> muss fuer ngen PTR allokiert sein, sowie stdgen[i] fuer 0 Zeilen;
   <stdtab> muss fuer ngen * 1 long integers mit stdtab[i][0] == 0.

   Koerperordnung und Zeilenlaenge muessen vor dem Aufruf gesetzt
   werden.

   Return: neue Dimension.
   ------------------------------------------------------------------ */
#define OPVEC(i) op_table[2*(i)]
#define OPGEN(i) op_table[2*(i)+1]


int zgensbasis(PTR seed, int noc, int seedcount, int ngen, Matrix_t *gen[], PTR space, 
    int *piv_table, PTR basis, int partdim, long *op_table, PTR stdgen[],
    long *std_tab[])

{
   static Matrix_t *transf = NULL;
   static long gencount = 1;
   long i, j, k;
   PTR xi, yi, xk, yk, temp, row, transfptr;
   FEL f;
   int igen;

   /* Initialize
      ---------- */
   if (transf == NULL)        /* identity matrix + one zero-row for the operation */
   {
      transf = matAlloc(ffOrder, noc + 1, noc);
      for (i = 0, row = transf->Data; i < noc; ++i)
      {
         ffInsert(row, i, FF_ONE);
         ffStepPtr(&row, noc);
      }
   }
   i = 1; 
   j = partdim + 1; 
   xi = ffGetPtr(space,partdim,noc);
   yi = ffGetPtr(basis,partdim,noc);
   k = partdim + 1; 
   xk = ffGetPtr(space,partdim,noc);
   yk = ffGetPtr(basis,partdim,noc);
   igen = 0;
   seed = ffGetPtr(seed,seedcount-1,noc);

   /* Main loop
      --------- */
   ffCopyRow(yk,seed, noc);
   ffCopyRow(xk,seed, noc);
   if (op_table != NULL)
   {
      OPVEC(k) = gencount;
      OPGEN(k) = 0;
   }
   myzcleanrow(xk, space, transf->Data, partdim, noc, piv_table);
   if ((piv_table[partdim] = ffFindPivot(xk, &f, noc)) < 0)
   {
      transfptr = matGetPtr(transf,partdim);
      ffMulRow(transfptr,FF_ZERO, noc);
      if (partdim < noc) 
         ffInsert(transfptr, partdim, FF_ONE);
      return partdim;
   } 
   gencount++;
   k++;
   partdim++;
   ffStepPtr(&xk, noc);
   ffStepPtr(&yk, noc);

   while (xi < xk)
   {
      ffMapRow(yi, gen[igen]->Data, noc, noc, yk);
      ffCopyRow(xk,yk, noc);

      /* Clean and check if we got a new vector
         -------------------------------------- */
      myzcleanrow(xk, space, transf->Data, partdim, noc, piv_table);
      if ((piv_table[partdim] = ffFindPivot(xk, &f, noc)) >= 0)
      {
         if (op_table != NULL)
         {
            OPVEC(k) = j;
            OPGEN(k) = igen+1;
         }
         ++k;
         partdim++;
         ffStepPtr(&xk, noc);
         ffStepPtr(&yk, noc);
      }
      else
      {
         std_tab[igen] = NREALLOC(std_tab[igen],long,++std_tab[igen][0] + 1);
         std_tab[igen][std_tab[igen][0]] = i;
         temp = ffAlloc(std_tab[igen][0], noc);
         memcpy(temp,stdgen[igen],ffRowSize(noc) * (std_tab[igen][0] - 1));
         row = ffGetPtr(temp,std_tab[igen][0] - 1, noc);
         transfptr = ffGetPtr(transf->Data,partdim, noc);
         ffCopyRow(row,transfptr, noc);
         if (partdim < noc) 
            ffInsert(row, partdim, FF_ZERO);
         ffMulRow(row, ffNeg(FF_ONE), noc);
         ffMulRow(transfptr, FF_ZERO, noc);
         if (partdim < noc) 
            ffInsert(transfptr, partdim, FF_ONE);
         sysFree(stdgen[igen]);
         stdgen[igen] = temp;
      }

      if (++igen >= ngen)        /* All the generators have been used */
      {
         igen = 0;
         ++i;
         ++j;
         ffStepPtr(&xi, noc);
         ffStepPtr(&yi, noc);
      }

   }
   return partdim;
}


long independent(Matrix_t *bas[], Matrix_t *mat, int dim, int **piv_table, 
    const long nummodgens, PTR dep)

{
    int i, j;
    FEL f;
    PTR matptr, basptr;
    MESSAGE(1,("independent: dim=%d\n",dim));
    for (i = 0; i < dim; i++)
    {
        if (bas[i] == NULL)
            continue;

        basptr = matGetPtr(bas[i],piv_table[i][0]);
        matptr = matGetPtr(mat,piv_table[i][0]);
        f = ffExtract(matptr, piv_table[i][1]);
        f = ffDiv(f, ffExtract(basptr, piv_table[i][1]));
        if (dep != NULL) {
            ffInsert(dep,i,f);
	}
        matAddMul(mat,bas[i],ffNeg(f));
    }
    piv_table[dim][0] = -1;
    if (nummodgens == -1 || big)
    {
        matptr = mat->Data;
        for (j = 0; j < mat->Nor && piv_table[dim][0] < 0; j++)
        {
            if ((piv_table[dim][1] = ffFindPivot(matptr, &f, mat->Noc)) >= 0)
                piv_table[dim][0] = j;
            ffStepPtr(&matptr,mat->Noc);
        }
    }
    else
    {
        int row = 0;
        matptr = mat->Data;
        for (j = 0; j < nummodgens && piv_table[dim][0] < 0; j++)
        {
            if ((piv_table[dim][1] = ffFindPivot(matptr, &f, mat->Noc)) >= 0)
                piv_table[dim][0] = row;
           matptr = ffGetPtr(matptr, dims[j], mat->Noc);
           row += dims[j];
        }
    }

    if (piv_table[dim][0] >= 0 && dep != NULL)
        ffInsert(dep, dim, FF_ONE);
    MESSAGE(2,("independent(): result=%d\n",piv_table[dim][0] >= 0));
    return piv_table[dim][0] >= 0;
}



Matrix_t *SmallForm(Matrix_t *mat)
{
    Matrix_t *small;
    int i, k;

    if ((small = matAlloc(mat->Field, numMgens + 1, dimM)) == NULL)
        return NULL;
    for(i = 0, k = 0; i <= numMgens; k += dims[i++])
        ffCopyRow(matGetPtr(small, i), matGetPtr(mat, k), dimM);
    matFree(mat);
    return small;
}


Matrix_t *bigform(Matrix_t *mat, Matrix_t **gens, long *op_table)
{
    long ind, max;
    Matrix_t *big;
    PTR bigptr, matptr, ptr;

    matptr = mat->Data;
    big = matAlloc(mat->Field, dimM, mat->Noc);
    bigptr = big->Data;
    max = 2 * dimM;
    for (ind = 2; ind <= max; ind += 2)
    {
        if (op_table[ind + 1] == 0)
        {
            ffCopyRow(bigptr, matptr, mat->Noc);
            ffStepPtr(&matptr, mat->Noc);
        }
        else
        {
            ptr = matGetPtr(big,op_table[ind] - 1);
            ffMapRow(ptr, gens[op_table[ind + 1] - 1]->Data, gens[0]->Nor, mat->Noc, bigptr);
        }
        ffStepPtr(&bigptr, mat->Noc);
    }
    return big;
}




Matrix_t **ringgens(Matrix_t *basis[], long n, long nummodgens, Matrix_t *regrep[], 
    char side, int big, Matrix_t **stdbas, long *op_table, Matrix_t **Ngen)

{
    int max = 0, next = 0, i, j, dim = 0, 
        **piv_table, **bpiv, **baspiv, d, g, *genind;
    Matrix_t **gens, *mat;
    PTR *regptr;
    FEL coeff;
    char name[100];

    if (side != 'l' && side != 'r')
    {
        mtxAbort(MTX_HERE,"Invalid side='%c'",side);
        return NULL;
    }
    if (   (piv_table = NALLOC(int *,n + 1)) == NULL
        || (baspiv = NALLOC(int *,n + 1)) == NULL)
        return NULL;

    for (i = 0; i <= n; i++)
    {
        if (   (piv_table[i] = NALLOC(int,2)) == NULL
            || (baspiv[i] = NALLOC(int,2)) == NULL)
            return NULL;
    }
    if (   (genind = NALLOC(int,n)) == NULL
        || (gens = NALLOC(Matrix_t *,n + 1)) == NULL
        || (regptr = NALLOC(PTR,n)) == NULL
        || (bpiv = NALLOC(int *,2)) == NULL
        || (bpiv[0] = NALLOC(int,2)) == NULL
        || (bpiv[1] = NALLOC(int,2)) == NULL)
        return NULL;

/* -----------------------------
   makes a basis for the algebra
   ----------------------------- */

    d = basis[0]->Noc;
    g = basis[0]->Nor;

    while (dim < n)
    {
            MESSAGE(1,("ringgens(): dim=%d\n",dim));
        if ((stdbas[dim] = matAlloc(ffOrder, g, d)) == NULL)
            return NULL;

/* chosing a random element of the algebra
   --------------------------------------- */
        for (i = 0; i < n; i++)
        {
            coeff = ffFromInt(mtxRandomInt(ffOrder));
            if (basis[i] != NULL)
                matAddMul(stdbas[dim], basis[i], coeff);
        }


/* testing if it is independent from the others
   -------------------------------------------- */

        if (!independent(stdbas, stdbas[dim], dim, piv_table, nummodgens, NULL))
        {
            matFree(stdbas[dim]);
            continue;
        }
        genind[max] = dim;
        if (big)
            gens[max] = bigform(stdbas[dim], Ngen, op_table);
        else    
            gens[max] = stdbas[dim];
        dim++;
	sprintf(name, "%s.%d", HomName, dim);
	matSave(gens[max], name);
        MESSAGE(1,("ringgens(): new element, dim=%d\n",dim));
        if ((regrep[max] = matAlloc(ffOrder, n, n)) == NULL)
            return NULL;
        regptr[max] = regrep[max]->Data;
        for (i = 0; i < genind[max]; i++)
        {
            if (side == 'r')
            {
                stdbas[dim] = matDup(stdbas[i]);
                matMul(stdbas[dim], gens[max]);
            }
            else         /* side == 'l' -- multiple from left */
            {
                stdbas[dim] = matDup(stdbas[genind[max]]);
                if (big)
                    mat = bigform(stdbas[i], Ngen, op_table);
                else
                    mat = matDup(stdbas[i]);
                matMul(stdbas[dim], mat);
                matFree(mat);
            }


            if (independent(stdbas, stdbas[dim], dim, piv_table, nummodgens, regptr[max]))
            {
		sprintf(name, "%s.%d", HomName, dim + 1);
		mat = (big == 0 ? matDup(stdbas[dim]) : bigform(stdbas[dim], Ngen, op_table));
		matSave(mat, name);
		matFree(mat);
                dim++;
                MESSAGE(1,("ringgens(): new element2, dim=%d\n",dim));
            }
            else
                matFree(stdbas[dim]);

            ffStepPtr(&(regptr[max]), n);
        }

        for (i = genind[max]; i < dim; i++)
        {
            Matrix_t *bigmat = NULL;
            if (side == 'l' && big)
                bigmat = bigform(stdbas[i], Ngen, op_table);
            for (next = 0; next <= max; next++)
            {
                if (side == 'r')
                {
                    stdbas[dim] = matDup(stdbas[i]);
                    matMul(stdbas[dim], gens[next]);
                }
                else /* side == 'l' -- multiple from left */
                {
                    stdbas[dim] = matDup(stdbas[genind[next]]);
                    if (big)
                        matMul(stdbas[dim], bigmat);
                    else
                        matMul(stdbas[dim], stdbas[i]);
                }

                if (independent(stdbas, stdbas[dim], dim, piv_table, nummodgens, regptr[next]))
                {
		    sprintf(name, "%s.%d", HomName, dim + 1);
		    mat = (big == 0 ? matDup(stdbas[dim]) : bigform(stdbas[dim], Ngen, op_table));
		    matSave(mat, name);
                    dim++;
                    MESSAGE(1,("ringgens(): new element3, dim=%d\n",dim));
                }
                else
                    matFree(stdbas[dim]);

                ffStepPtr(&(regptr[next]),n);
            }
            if (big && side == 'l')
                matFree(bigmat);
        }
        max++;
    }

    if (!big)
    {
        for (j = 0; j < n; j++)
            matFree(basis[j]);
    }

    if (side == 'l')        /* left repr. must be transposed */
    {
        for (i = 0; i < max; i++)
        {
            mat = matTransposed(regrep[i]);
            matFree(regrep[i]);
            regrep[i] = mat;
        }
    }
    regrep[max] = NULL;
    gens[max] = NULL;
    return (gens);
}



/* ========================================================================== */






/* ========================================================================== */





static int ParseArgs()

{
    int tmp;

    /* Options.
       -------- */
    standard = appGetOption(App,"-t");
    hominstd = appGetOption(App,"-s");
    tmp = appGetOption(App,"-r");
    if (tmp)
	reg = 'r';
    if (appGetOption(App, "-l"))
	reg = 'l';
    if (tmp && reg == 'l')
    {
	mtxAbort(MTX_HERE,"-r and -l cannot be used simultaneously");
	return -1;
    }
    if (reg == 'r' || reg == 'l')
    {
        hominstd = 1;
        standard = 1;
    }
    big = appGetOption(App,"-b");
    if(big != 0 && reg == '?')
        {
                MESSAGE(0, ("-b is only used with -r/l\n"));
                big = 0;
        }

    hd = appGetIntOption(App,"-H",0,1,1000000);


    /* Arguments.
       ---------- */
    if (appGetArguments(App,3,3) < 0)
        return -1;
    MName = App->ArgV[0];
    NName = App->ArgV[1];
    HomName = App->ArgV[2];
    comp = strcmp(MName,NName);
    if (hominstd && comp)
    {
        mtxAbort(MTX_HERE,"-s requires <M> = <N>");
        return -1;
    }

    return 0;
}



static int ReadFiles()
{
    /* Read the .cfinfo files.
       ----------------------- */
    if (latReadInfo(&MInfo,MName) != 0)
        return -1;

    /* Load the generators.
       -------------------- */
    MESSAGE(1,("Reading generators\n"));
    if ((MRep = mrLoad(MInfo.BaseName,MInfo.NGen)) == NULL)
        return -1;
    dimM = MRep->Gen[0]->Noc;
    if (comp)
    {
        if ((NRep = mrLoad(NName,MInfo.NGen)) == NULL)
            return -1;
    }
    else
        NRep = MRep;

    /* Read the head, if -H is used.
       ----------------------------- */
    if (hd > 0)
    {
        Matrix_t *tmp;
        char fn[200];
        sprintf(fn, "%s.rad",MName);
        MESSAGE(1,("Reading the head (%s)\n",fn));
        if((tmp = matLoad(fn)) == NULL)
            return -1;
        if((rad = matCutRows(tmp, hd, dimM - hd)) == NULL)
            return -1;
        matEchelonize(rad);
        rad->Data = (PTR) sysRealloc(rad->Data, ffSize(dimM, tmp->Noc));
        matFree(tmp);
    }

    return 0;
}





static int AllocateWorkspace()
{
    int i;

    if (   (basis = ffAlloc(dimM + 1, dimM)) == NULL
        || (space = ffAlloc(dimM + 1, dimM)) == NULL)
        return -1;
    if ((piv = NALLOC(int,dimM + 2)) == NULL)
        return -1;
    if ((op = NALLOC(long,2 * dimM + 2)) == NULL)
        return -1;
    if (   (stdgen = NALLOC(PTR,MInfo.NGen)) == NULL 
        || (stdtab = NALLOC(long *,MInfo.NGen)) == NULL)
        return -1;
    for (i = 0; i < MInfo.NGen; i++)
    {
        if (   (stdgen[i] = ffAlloc(0, dimM)) == NULL
            || (stdtab[i] = ALLOC(long)) == NULL)
            return -1;
        stdtab[i][0] = 0;
    }

/* NEW PART
   -------- */
    for (i = 0; i < MInfo.NCf; i++)
        maxnumMgens += MInfo.Cf[i].mult;
    if ((posimages = NALLOC(Matrix_t **, maxnumMgens)) == NULL
    || (dims = NALLOC(int, maxnumMgens)) == NULL
    || (kerdim= NALLOC(int, maxnumMgens)) == NULL
    || (esyspiv = NALLOC(int, maxnumMgens* NRep->Gen[0]->Nor)) == NULL
    || (esys = matAlloc(ffOrder, 0, 0)) == NULL)
        return 1;

/* END OF NEW PART
   --------------- */
    return 0;
}


static int Init(int argc, char **argv)

{
    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
        return -1;
    if (ParseArgs() != 0)
        return -1;
    if (ReadFiles() != 0)
        return -1;
    if (AllocateWorkspace() != 0)
    {
        mtxAbort(MTX_HERE,"Cannot allocate work space");
        return -1;
    }
    return 0;
}


static void Cleanup()

{
    if (MRep != NULL)
        mrFree(MRep);
    if (NRep != NULL && NRep != MRep)
        mrFree(NRep);
    appFree(App);
}







/* =========================================================================
   spinpartstdbas() spin up a newdim dimensional part of the 
   spinning basis generated by vec beginning at partdim
   ========================================================================= */

Matrix_t *spinpartstdbas(PTR vec, const long *op_table, Matrix_t *gens[], 
    int part_dim, int newdim)

{
   MTX_ASSERT(gens != NULL && gens[0] != NULL, NULL);
    Matrix_t *mat;
    PTR ptr, row;
    int l, newpartdim;

    newpartdim = newdim + part_dim;

    mat = matAlloc(ffOrder,newdim,gens[0]->Noc);
    ptr = mat->Data;
    ffCopyRow(ptr,vec, gens[0]->Noc);
    ffStepPtr(&ptr, gens[0]->Noc);
    for (l = part_dim + 2; l <= newpartdim; l++)
    {
        row = matGetPtr(mat,op_table[2*l] - 1 - part_dim);
        ffMapRow(row, gens[op_table[2*l + 1] - 1]->Data, gens[0]->Nor, gens[0]->Noc,ptr);
        ffStepPtr(&ptr, gens[0]->Noc);
    }
    return mat;
}

/* =========================================================================
   veccont() checks if vec is contained in the subspace generated by mat
   ========================================================================= */

int veccont(Matrix_t *mat, PTR vec, int *pivot_table)
{
    PTR v;
    FEL f;
    int is_contained;

    v = ffAlloc(1, mat->Noc);
    ffCopyRow(v,vec, mat->Noc);
    ffCleanRow(v,mat->Data,mat->Nor,mat->Noc,pivot_table);
    is_contained = ffFindPivot(v,&f, mat->Noc) < 0;
    sysFree(v);
    return is_contained;
}




static int MakeKernels(int cf, Matrix_t **ker1, Matrix_t **ker2)

{
    char file_name[200];
    long t0 = sysTimeUsed();

    /* Load the peak word kernel for M.
       -------------------------------- */
    sprintf(file_name, "%s%s.k", MName, latCfName(&MInfo,cf));
    if ((ker1 != NULL) && ((*ker1 = matLoad(file_name)) == NULL))
    {
        mtxAbort(MTX_HERE,"Cannot load %s -- did you run 'pwkond %s'?",
            file_name,MName);
        return 1;
    }

    if (ker2 == NULL)
        return 0;
    /* If N is different from M, find the stable peak word kernel in N.
       ---------------------------------------------------------------- */
    if (!comp)
        *ker2 = *ker1;
    else
    {
        Matrix_t *word2;
        WgData_t *wg = wgAlloc(NRep);
        MESSAGE(1,("Calculating the stable peak word kernel in %s\n",NName));
        word2 = wgMakeWord(wg,MInfo.Cf[cf].peakword);
        wgFree(wg);
        matInsert_(word2, MInfo.Cf[cf].peakpol);

        StablePower_(word2,NULL,ker2);
        matFree(word2);
    }
    MESSAGE(0, ("                %ld\n", sysTimeUsed() - t0));
    return 0;
}


/* ------------------------------------------------------------------------ */

int main(int argc, char **argv)

{
    int rc = 0;                        /* Program exit code */
    long t0, tposim = 0, teqs = 0, tstker = 0, tgauss = 0,
                    tspbas = 0;        /* Time */
    Matrix_t 
             *spinbas,                /* the spinning basis */
             *spinbasi = 0,         /* its inverse) */
             *ker1,                /* nullspace in M */
             *ker2 = NULL,        /* nullspace in N */
             *mat,
             *result,
             *partesys,
             *tresys,
             **gens,
             **stdbas = NULL,
             *echker = NULL,
             **homs,
             **regrep;
    PTR basptr, kerptr, echkerptr,
        stdgenptr, sysptr, resptr, partptr, esysptr;
    int sb, seedcount = 0, row, oldNor = 0, m, s, k,
        newpartdim, ind, col;
    char name[100];
    int i, l;
    FEL f;
    Lat_Info end_info;


    if (Init(argc,argv) != 0)
    {
        mtxAbort(MTX_HERE,"Initialization failed");
        return 1;
    }


    /* Main loop: for each constituent of M.
       ------------------------------------- */
    for (i = 0; i < MInfo.NCf; i++)
    {
        int j;

        MESSAGE(0,("Next constituent: %s%s\n",MName, latCfName(&MInfo,i)));

        t0 = sysTimeUsed();
        if (comp && (MakeKernels(i,&ker1,NULL) != 0))  /* Kernels of the peakwords */
            return 1;
        if (!comp && MakeKernels(i,&ker1,&ker2) != 0)  /* Kernels of the peakwords */
            return 1;
        tstker += sysTimeUsed() - t0;
        seedcount = 0;
        if (hd)
        {
            echker = matDup(ker1);
            echkerptr = echker->Data;
        }

        /* Make the next part of the spinning basis in M.
           ---------------------------------------------- */
        for (j = 0; j < ker1->Nor; j++)
        {
            t0 = sysTimeUsed();
            seedcount++;
            MESSAGE(1,("Taking kernel vector %d\n",j+1));
            if (hd)
            {
                ffCleanRow(echkerptr, rad->Data, rad->Nor, dimM, rad->PivotTable);
                if((ffFindPivot(echkerptr, &f, dimM)) < 0)
                {
                    ffStepPtr(&echkerptr, dimM);
                    continue;
                }
                ffStepPtr(&echkerptr, dimM);
            }
            if ((newpartdim = zgensbasis(ker1->Data, dimM, seedcount, MInfo.NGen, MRep->Gen, 
                space, piv, basis, partdim, op, stdgen, stdtab)) == partdim) 
            {
                MESSAGE(1,("No new basis vectors - skipping\n"));
                continue;
            }
            MESSAGE(0,("Vector %d (seedcount=%d) spins up to %d\n",j+1,
                seedcount,newpartdim));
            dims[numMgens] = newpartdim - partdim;
            tspbas += sysTimeUsed() - t0;
           
            /* Extending the radical with the new part of the module */

            if(hd > 0 && newpartdim < dimM)
            {
                PTR radptr;
                basptr = ffGetPtr(basis, partdim, dimM);
                radptr = ffGetPtr(rad->Data, rad->Nor, dimM);
                for(k = 0; k < dims[numMgens]; k++)
                {
                    ffCopyRow(radptr, basptr, dimM);
                    ffCleanRow(radptr, rad->Data, rad->Nor, dimM, rad->PivotTable);
                    rad->PivotTable[rad->Nor] = ffFindPivot(radptr, &f, dimM);
                    if(rad->PivotTable[rad->Nor] >= 0)
                    {
                        rad->Nor++;
                        ffStepPtr(&radptr, dimM);
                    }
                    ffStepPtr(&basptr, dimM);
                }
            }

            t0 = sysTimeUsed();
            if (ker2 == NULL)
                MakeKernels(i, NULL, &ker2);
            tstker += sysTimeUsed() - t0;
            kerdim[numMgens] = ker2->Nor;

/* -----------------------------------------------
   making the possible images in the second module
   ----------------------------------------------- */

            t0 = sysTimeUsed();


/* MODIFIED PART
   ------------- */

            MESSAGE(1,("Calculating the possible images in %s\n",NName));
            if ((posimages[numMgens] = NALLOC(Matrix_t *, ker2->Nor)) == NULL)
                return 1;
            kerptr = ker2->Data;
            for (k = 0; k < ker2->Nor; k++)
            {
                posimages[numMgens][k] =
                   spinpartstdbas(kerptr, op, NRep->Gen, partdim, dims[numMgens]);
                ffStepPtr(&kerptr, ker2->Noc);
            }

/* END OF MODIFIED PART
   -------------------- */
            
            tposim += sysTimeUsed() - t0;


/* ------------------------------------------------------------------
   builds up the part of the system of equations corresponding to the
   relations obtained at the last step.
   ------------------------------------------------------------------ */

            if ((mat = matAlloc(ffOrder, esys->Nor + ker2->Nor, esys->Noc + ker2->Nor)) == NULL)
                return 1;
            matMulScalar(mat, FF_ZERO);
            if (matCopyRegion(mat, 0, 0, esys, 0, 0, oldNor, esys->Noc) == -1)
                return 1;
            matFree(esys);
            esys = mat;
            MESSAGE(1,("Building equation system (%dx%d)\n",
                    NRep->Gen[0]->Noc * dims[numMgens] * MInfo.NGen, mat->Nor));

            if (esys->Nor == 0)   /* there are no homomorphisms */
            {
                if (newpartdim == MRep->Gen[0]->Nor)
                {
                    MESSAGE(0,("Warning: There are no homomorphisms from "
                        "%s to %s\n", MName, NName));
                   return 0;
                }
                partdim = newpartdim;
                for (k = 0; k < MInfo.NGen; k++)
                {
                    sysFree(stdgen[k]);
                    if ((stdgen[k] = ffAlloc(0, esys->Noc)) == NULL)
                        return 1;
                    stdtab[k][0] = 0;
                }
		numMgens++;
                continue;
            }
            esysptr = ffGetPtr(esys->Data, oldNor, esys->Noc);


/* builds up the system of equations
   --------------------------------- */

            if ((tresys = matAlloc(ffOrder, esys->Noc, NRep->Gen[0]->Nor)) == NULL)
                return 1;
            for (k = 0; k < MInfo.NGen; k++)        /* loop for the generators */
            {
                stdgenptr = stdgen[k];

                for (l = 1; l <= stdtab[k][0]; l++)
                                        /* loop for the vectors */
                {
                    int t;
                    t0 = sysTimeUsed();
                    matMulScalar(tresys, FF_ZERO);
                    // the equations for one vector
                    sysptr = tresys->Data;
                    col = 0;
                    for (m = 0; m <= numMgens; m++)
                    {
                        for (s = 0; s < kerdim[m]; s++)
                        {
                            if (m == numMgens)
                            {
                                basptr = matGetPtr(posimages[m][s], stdtab[k][l] - 1);
                                ffMapRow(basptr, NRep->Gen[k]->Data, NRep->Gen[0]->Nor, NRep->Gen[0]->Noc, sysptr);
                                ffMulRow(sysptr, ffNeg(FF_ONE), NRep->Gen[0]->Noc);
                            }
                            basptr = posimages[m][s]->Data;
                            for (sb = 0; sb < dims[m]; sb++)
                            {
                                if ((f = ffExtract(stdgenptr, col + sb)) != FF_ZERO)
                                    ffAddMulRow(sysptr, basptr, f, NRep->Gen[0]->Noc);
                                ffStepPtr(&basptr, NRep->Gen[0]->Noc);
                            }
                            ffStepPtr(&sysptr, NRep->Gen[0]->Noc);
                        } 
                        col += dims[m];
                    }
                    ffStepPtr(&stdgenptr, dimM);



                    teqs += sysTimeUsed() - t0;

                    // eliminating the superfluous equations

                    t0 = sysTimeUsed();
                    partesys = matTransposed(tresys);
                    partptr = partesys->Data;
                    for (t = 0; t < partesys->Nor; t++)
                    {
                        ffCleanRow(partptr, esys->Data , oldNor, partesys->Noc, esyspiv);

                        if ((esyspiv[oldNor] = ffFindPivot(partptr, &f, partesys->Noc)) >= 0)
                        {
                            ffCopyRow(esysptr, partptr, partesys->Noc);
                            if ((++oldNor) > esys->Noc) 
                            {
                                mtxAbort(MTX_HERE,"The matrix has rank greater than number of rows");
                                return 1;
                            }
                            ffStepPtr(&esysptr, partesys->Noc);
                        }
                        ffStepPtr(&partptr, partesys->Noc);
                    }
                    matFree(partesys);
                    tgauss += sysTimeUsed() - t0;
                }        /* end of loop with l */
            }                /* end of loop with k */

            matFree(tresys);

            MESSAGE(1, ("%d homomorphisms found\n", esys->Noc - oldNor));

/* gives back the superfluous space to the memory 
   ---------------------------------------------- */
            partdim = newpartdim;
            for (k = 0; k < MInfo.NGen; k++)
            {
                sysFree(stdgen[k]);
                if ((stdgen[k] = ffAlloc(0, esys->Noc)) == NULL)
                    return 1;
                stdtab[k][0] = 0;
            }


            if (newpartdim == dimM)
            {
                spinbas = matAlloc(ffOrder,dimM,dimM);
                sysFree(spinbas->Data);
                spinbas->Data = basis;
                sprintf(name,"%s.spb",MName);
                MESSAGE(1, ("Writing spinning basis to %s\n", name));
                matSave(spinbas, name);
                if (standard || hominstd)
                    spinbasi = matInverse(spinbas);
                if (standard)
                {
                    MESSAGE(1,("Transforming %s into spinning basis\n",
                        MName));
                    for (k = 0; k < MInfo.NGen; k++)
                    {
                        mat = matDup(spinbas);
                        matMul(mat, MRep->Gen[k]);
                        matMul(mat, spinbasi);
                        sprintf(name, "%s.std.%d", MName, k + 1);
                        matSave(mat, name);
                        if (reg != '?')
                        {
                            matFree(NRep->Gen[k]);
                            NRep->Gen[k] = mat;
                        }
                        if (reg == '?')
                            matFree(mat);
                    }
                }
/* solving the system of equations */


                if ((tresys = matTransposed(esys)) == NULL)
                    return 1;
                if ((result = matNullSpace_(tresys,0)) == NULL)
                    return 1;
                homs = NALLOC(Matrix_t *, result->Nor);
                resptr = result->Data;
                for (row = 0; row < result->Nor; row++)
                {
                    if ((homs[row] = matAlloc(ffOrder, dimM, NRep->Gen[0]->Nor)) == NULL)
                        return 1;
                    col = 0;
                    k = 0;
                    for (m = 0; m <= numMgens; m++) {
                        if ((mat = matAlloc(ffOrder, dims[m], NRep->Gen[0]->Noc)) == NULL)
                            return 1;
                        matMulScalar(mat, FF_ZERO);
                        for (ind = 0; ind < kerdim[m]; ind++)
                        {
                            if ((f = ffExtract(resptr, col + ind)) != FF_ZERO) {
                                matAddMul(mat, posimages[m][ind], f);
                            }
                        }
                        matCopyRegion(homs[row], k, 0, mat, 0, 0, -1, -1);
                        col += kerdim[m];
                        k += dims[m];
                        matFree(mat);
                    }
                    if (hominstd)
                        matMul(homs[row], spinbasi);
                    if(reg == '?')
		    {
                        sprintf(name, "%s.%d", HomName, row + 1);
                        matSave(homs[row], name);
                        matFree(homs[row]);
		    }
                    if(big)
                        homs[row] = SmallForm(homs[row]);
                    ffStepPtr(&resptr, result->Noc);
                }
                if (reg == '?')
                    return rc;

                MESSAGE(1,("Calculating regular representation\n"));
                if ((regrep = NALLOC(Matrix_t *, result->Nor)) == NULL)
                    return 1;
                if ((stdbas = NALLOC(Matrix_t *, result->Nor + 1)) == NULL)
                    return 1;

                gens = ringgens(homs, result->Nor, numMgens + 1, regrep, reg, big, stdbas, op, NRep->Gen);

                for (k = 0; gens[k] != NULL; k++)
                {
                            sprintf(name, "%s.gens.%d", HomName, k + 1);
                    matSave(gens[k], name);
                    sprintf(name, "%s.%crr.%d", HomName, reg, k + 1);
                    matSave(regrep[k], name);
                }

                /* Cretae the <endo>.lrr.cfinfo file */
                memset(&end_info,0,sizeof(end_info));
                end_info.NGen = k;
                sprintf(end_info.BaseName,"%s.%crr",HomName,reg);
                latWriteInfo(&end_info);
                return rc;
            }                /* endif dimM == newpartdim */
            numMgens++;
        }                /* end of loop with j */
        if (ker2 != NULL)
        {
            matFree(ker2);
            ker2 = NULL;
        }
        if (comp) 
            matFree(ker1);
    }


    Cleanup();
    return rc;
}


/**
@page prog_mkhom mkhom - Homomorphisms

@section mkhom_syntax Command Line
<pre>
mkhom [@em Options] [-ts] [-r|-l] [-b @em Mode] [-H @em Dim] @em M @em N @em Hom
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -t
  Calculate generators for @em M in spinning basis.
@par -s
  When @em M equals @em N, calculate endomorphisms in spinning basis.
@par -r
  When @em M equals @em N, find a generating set of End(M), and
  calculate the right regular representation.
@par -l
    Like -r but calculate the left regular representation.
@par -b @em Mode
    Save memory. @em Mode can be 0, 1, or 2 (see description).
@par -H @em Dim
    If the radical is given, @em Dim is the dimension of the head.
@par @em M
    First representation.
@par @em N
    Second representation.
@par @em Hom
    Homomorphisms.
  

@section mkhom_inp Input Files

@par @em M.cfinfo, @em N.cfinfo
  Constituent information.
@par @em M.1, @em M.2,...
  Generators in the first representation.
@par @em N.1, @em N.2,...
  Generators in the second representation.
@par @em M.rad
  Generators for the head of @em M (with -H).
@par @em MCf.k
  Uncondense matrix produced by @ref prog_pwkond "pwkond".

@section mkhom_out Output Files
@par @em M.std
  The spinning basis for @em M.
@par @em Hom.1, @em Hom.2,...
  A k-basis of Hom(@em M,@em N).
@par @em M.std.1, @em M.std.2,...
  Generators in spinning basis (with -t).

@section mkhom_desc Description
This program calculates a basis for the vector space of homomorphisms
between two kG-modules, Hom<sub>kG</sub>(M,N). In the case M=N the program
optionally finds a generating set for the algebra of endomoprhisms,
End<sub>kG</sub>(M), and calculates the corresponding left or right regular 
representation. 

If used without any options, @p mkhom writes the spinning basis of M
to @em M.std, and a k-basis of the homomorphism space to @em Hom.1, @em Hom.2,...
The latter are given with respect to the spinning basis of M and the original basis of N.
To get the homomorphisms with respect to the original bases of M
and N, multiply the matrices from the left with the inverse of @em Hom.

@p mkhom uses peak words of the first module. Thus, before using the program, 
@ref prog_chop "chop" and @ref prog_pwkond "pwkond" must have been run on the first module.

@section mkhom_impl Implementation Details
The algorithm used by this program was developed by Magdolna Szöke [@ref Sz98].
**/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
