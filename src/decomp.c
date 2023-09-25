////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Decompose a module into direct summands.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

static MtxApplicationInfo_t AppInfo = {
"decomp", "Decompose module",
"SYNTAX\n"
"    decomp [-QVta] <M> <Endo>\n"
"\n"
"ARGUMENTS\n"
"    <Module> ................ Module to decompose.\n"
"    <Endo> .................. Endomorphism ring.\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -t ...................... Write transformed generators.\n"
"    -a ...................... Write the action on the direct summands.\n"
"\n"
"FILES\n"
"    <M>.{1,2...} ............ I  Generators on <M>.\n"
"    <M>.cfinfo .............. I  Constituent info file for <M>.\n"
"    <Endo>.{1,2...} ......... I  k-Basis of the endomorphism ring.\n"
"    <Endo>.gens.{1,2...} .... I  Generating system of the endomorphism ring.\n"
"    <Endo>.lrr.{1,2...} ..... I  Left regular repr. of the endomorphism ring.\n"
"    <Endo>.lrr.cfinfo ....... I  Constituent info file for <Endo>.lrr after\n"
"                                 running CHOP and PWKOND.\n"
"    <Endo>.lrr.soc .......... I  Basis of the socle of <Endo>.lrr (made by SOC)\n"
"    <M>.dec ................. O  Basis of <M> reflecting the decomposition.\n"
"    <M>.dec.{1,2...} ........ O  Generators in decomp. basis (with -t).\n"
"    <M>.<Comp>.{1,2...} ..... O  Generators on the components (with -a).\n"
};                                                                              
                                                                                
static MtxApplication_t *App = NULL;                                            
static const char *ModName = NULL;
static const char *EndoName = NULL;
static Lat_Info ModInfo;		/* Data from .cfinfo file */
static Lat_Info LrrInfo;		/* Data from .cfinfo file */
static int moddim = 0;
static int enddim = 0;
static int headdim = 0;
static int compdim[LAT_MAXCF];
static char compnm[LAT_MAXCF];
static Matrix_t *head = NULL;
static int TransformGenerators = 0;	/* -t: Transform into decomp. basis */
static int WriteAction = 0;		/* -a: Write action on components */

////////////////////////////////////////////////////////////////////////////////////////////////////

static void ParseArgs()
{
    TransformGenerators = appGetOption(App,"-t");
    WriteAction = appGetOption(App,"-a");
    appGetArguments(App,2,2);
    ModName = App->argV[0];
    EndoName = App->argV[1];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void ReadFiles()
{
    char fn[200];
    Matrix_t *tmp, *tmp2;
    int i;

    /* Read the .cfinfo files and calculate some dimensions.
       ----------------------------------------------------- */
    latReadInfo(&ModInfo,ModName);
    sprintf(fn,"%s.lrr",EndoName);
    latReadInfo(&LrrInfo,fn);
    moddim = 0;
    for (i = 0; i < ModInfo.nCf; ++i)
	moddim += ModInfo.Cf[i].dim * ModInfo.Cf[i].mult;

    enddim = headdim = 0;
    for (i = 0; i < LrrInfo.nCf; i++)
    {
	enddim += LrrInfo.Cf[i].dim * LrrInfo.Cf[i].mult;
	headdim += LrrInfo.Cf[i].dim * LrrInfo.Cf[i].dim / LrrInfo.Cf[i].spl;
    }
    if (headdim > enddim || headdim <= 0)
    {
	mtxAbort(MTX_HERE,"The head (%d) is bigger than the ring itself (%d)!",
		headdim,enddim);
    }
    MESSAGE(1,("dim(M)=%d, dim(E)=%d, dim(Head)=%d\n",moddim,enddim,headdim));

    // Read the basis of the head.
    sprintf(fn,"%s.lrr.soc",EndoName);
    MESSAGE(1,("Loading socle basis\n"));
    tmp = matLoad(fn);
    tmp2 = matInverse(tmp);
    matFree(tmp);
    tmp = matTransposed(tmp2);
    matFree(tmp2);
    head = matCutRows(tmp,0,headdim);
}




static void Init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    ParseArgs();
    ReadFiles();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void WriteOutput(Matrix_t *bas)
{
    char name[200];
    MatRep_t *rep = NULL;
    int i;

    /* Write the decomposition basis.
       ------------------------------ */
    sprintf(name,"%s.dec", ModName);
    MESSAGE(1,("Writing the decomposition basis (%s)\n",name));
    matSave(bas,name);

    /* Transform the generators.
       ------------------------- */
    if (TransformGenerators || WriteAction)
    {
	MESSAGE(1,("Transforming the generators\n"));
	sprintf(name,"%s.std",ModName);
	rep = mrLoad(name, ModInfo.NGen);
	mrChangeBasis(rep,bas);
	if (TransformGenerators)
	{
	    sprintf(name,"%s.dec",ModName);
	    MESSAGE(1,("Writing transformed generators (%s.1, ...)\n",name));
	    mrSave(rep,name);
	}
    }


    /* Write the action of the generators on the direct summands.
       ---------------------------------------------------------- */
    if (WriteAction)
    {
	MESSAGE(1,("Writing the action on the direct summands\n"));
	for (i = 0; i < rep->NGen; i++)
	{
	    int k, block_start = 0;
	    for (k = 0; k < LrrInfo.nCf; k++)
	    {
		int l;
		for (l = 0; l < LrrInfo.Cf[k].dim / LrrInfo.Cf[k].spl; l++)
		{
		    Matrix_t *tmp = matCut(rep->Gen[i],block_start,block_start,
			compdim[k],compdim[k]);
		    block_start += compdim[k];
		    sprintf(name, "%s.comp%d%c%d.%d", ModName,compdim[k], 
			compnm[k],l+1,i+1);
		    matSave(tmp, name);
		    matFree(tmp);
		}
	    }
	}
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    int i, j, l, num, dim = 0;
    Matrix_t *bas, *partbas = NULL, *mat, *ker;
    PTR headptr;
    FEL f;
    char name[200];
    FPoly_t *pol = NULL;

    Init(argc,argv);

/* -------------------------------------------------------
   makes the corresponding element of the endomorphismring
   ------------------------------------------------------- */

    bas = matAlloc(ffOrder,moddim,moddim);
    headptr = head->data;
    for (i = 0; i < LrrInfo.nCf; i++)
    {
	MESSAGE(1,("Next constituent is %s%s\n",LrrInfo.BaseName,
	    latCfName(&LrrInfo,i)));
	for (j = 0; j < LrrInfo.Cf[i].dim / LrrInfo.Cf[i].spl; j++)
	{
	    num = LrrInfo.Cf[i].dim;
	    do
	    {
		if (partbas != NULL)
		    matFree(partbas);
		partbas = matAlloc(ffOrder, moddim, moddim);
		if (num-- == 0)
		{
		    mtxAbort(MTX_HERE,"na, most mi van?");
		    return 1;
		}
		for (l = 0; l < enddim; l++)
		{
		    if ((f = ffExtract(headptr, l)) == FF_ZERO)
		    	continue;
		    sprintf(name, "%s.%d", EndoName, l+1);
		    if ((mat = matLoad(name)) == NULL)
			return 1;
		    matAddMul(partbas,mat,f);
		    matFree(mat);
		}
		ffStepPtr(&headptr, enddim);

		if (pol != NULL)
		    fpFree(pol);
		pol = charpol(partbas);
	    }
	    while (LrrInfo.Cf[i].dim != 1 && pol->NFactors == 1 
		&& pol->Factor[0]->degree == 1
		&& pol->Factor[0]->data[0] == FF_ZERO 
		&& pol->Factor[0]->data[1] == FF_ONE); /* i.e.,charpol == x^enddim */
	    headptr = ffGetPtr(headptr,num, enddim);


            /* Make the stable kernel.
               ----------------------- */
	    StablePower_(partbas,NULL,&ker);
	    compdim[i] = moddim - ker->nor;
	    for (l = i - 1; l >= 0 && compdim[i] != compdim[l]; l--)
		;
	    if (l >= 0)
		compnm[i] = (char)(compnm[l] + 1);
	    else
		compnm[i] = 'a';
	    matFree(ker);
	    MESSAGE(0,("The %d-th direct summand is: %d%c\n\n", j,
		compdim[i], compnm[i]));


	    /* Append <partbas> to <bas>.
	       -------------------------- */
	    matEchelonize(partbas);
	    #ifdef MTX_DEBUG
	    #endif
	    matCopyRegion(bas,dim,0,partbas,0,0,-1,-1);
	    dim += partbas->nor;
	}
    }
    matFree(partbas);

    if (dim != moddim)
    {
	mtxAbort(MTX_HERE,"Something is wrong - dimension mismatch (%d vs. %d)", dim, moddim);
	for (i = 0; i < LrrInfo.nCf; i++)
    	    printf("%d  ", compdim[i]);
    	printf("\n");
	return 1;
    }



    WriteOutput(bas);


    appFree(App);
    return 0;
}



/**
@page prog_decomp decomp - Decompose a Module                                                  

@section decomp_syntax Command Line
<pre>
decomp @em Options [-ta] @em Module @em Endo @em RadBasis
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -t
    Transform generators.
@par -a
    Calculate the action on direct summands.
@par @em Module
    Name of the module to decompose.
@par @em Endo
    Name of the endomorphism ring.
@par @em RadBasis
    Basis for radical series of the endomorphism ring.

@section decomp_inp Input Files
@par @em Module.1, @em Module.2, ...
  Generators for the module.
@par @em Module.cfinfo
  Constituent information.
@par @em Endo.1, @em Endo.2, ...
  A k-Basis of the endomorphism ring.
@par @em Endo.gens.1, @em Endo.gens.2, ...
  A generating system of the endomorphism ring.
@par @em Endo.lrr.1, @em Endo.lrr.2, ...
  Left regular representation of the endomorphism ring.
@par @em Endo.lrr.cfinfo
  Constituent information for the left regular representation.
  At least @ref prog_chop "chop" and @ref prog_mksub "mksub"
  must have been run.
@par @em Endo.lrr.soc
  Basis if the socle of @em Endo.lrr (make by @ref prog_soc "soc").

@section decomp_out Output Files
@par @em Module.dec
  Basis of the module reflecting the decomposition.
@par @em Module.dec.1, @em Module.dec.2, ...
  Basis of the module reflecting the decomposition.
  Generators in decomp. basis (with -t).
@par @em Module.@em Comp.1, @em Module.@em Comp.2, ...
  Generators on the components (with -a).

@section decomp_desc Description
This program decomposes a module into its direct summands, using the
head of the endomorphism ring. It assumes that the endomorphism ring
and its left regular representation has already be calculated using
@ref prog_mkhom "mkhom". It also expects that @ref prog_chop "chop".
@ref prog_pwkond "pwkond", and @ref prog_rad "rad" have been run on
the left regular representation.

@b decomp produces three types of output files:
- A basis of the module reflecting the decomposition is written to 
  @em Module.dec. With respect to this basis, the generators have
  a block-diagonal structure with the blocks corresponding to the
  direct summands.
- If you use the "-t" option, @b decomp calculates the action of 
  the generators with respect to the decomposition
  basis, and writes it to @em Module.dec.1, @em Module.dec.2,...
- If you use the "-a" option, the program also calculates the action 
  of the generators on each direct summand, i.e., the blocks of the 
  matrices above. They are written to @em Module.@em Comp.1, @em Module.@em Comp.2,
  and so on, where @em Comp is the name 
  of the direct summand, containing the isomorphism type (dimension 
  plus one letter) and a number counting isomorphic summands that 
  occur more than once in the decomposition.
                                                                                
@section decomp_impl Implementation Details
  The algorithm used by this program was developed by Magdolna Szöke [@ref Sz98].

**/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
