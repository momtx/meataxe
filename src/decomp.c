/* ============================= C MeatAxe ==================================
   File:        $Id: decomp.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Decompose a module into direct summands.
   --------------------------------------------------------------------------
   Written by Magdolna Szoke.
   Revised by Michael Ringe.
   (C) Copyright 1999 Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"


MTX_DEFINE_FILE_INFO

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

static int ParseArgs()

{
    TransformGenerators = AppGetOption(App,"-t");
    WriteAction = AppGetOption(App,"-a");
    if (AppGetArguments(App,2,2) < 0)
	return -1;
    ModName = App->ArgV[0];
    EndoName = App->ArgV[1];
    return 0;
}


static int ReadFiles()

{
    char fn[200];
    Matrix_t *tmp, *tmp2;
    int i;

    /* Read the .cfinfo files and calculate some dimensions.
       ----------------------------------------------------- */
    if (Lat_ReadInfo(&ModInfo,ModName) != 0)
        return -1;
    sprintf(fn,"%s.lrr",EndoName);
    if (Lat_ReadInfo(&LrrInfo,fn) != 0)
        return -1;
    moddim = 0;
    for (i = 0; i < ModInfo.NCf; ++i)
	moddim += ModInfo.Cf[i].dim * ModInfo.Cf[i].mult;

    enddim = headdim = 0;
    for (i = 0; i < LrrInfo.NCf; i++)
    {
	enddim += LrrInfo.Cf[i].dim * LrrInfo.Cf[i].mult;
	headdim += LrrInfo.Cf[i].dim * LrrInfo.Cf[i].dim / LrrInfo.Cf[i].spl;
    }
    if (headdim > enddim || headdim <= 0)
    {
	MTX_ERROR2("The head (%d) is bigger than the ring itself (%d)!",
		headdim,enddim);
	return -1;
    }
    MESSAGE(1,("dim(M)=%d, dim(E)=%d, dim(Head)=%d\n",moddim,enddim,headdim));

    /* Read the basis of the head.
       --------------------------- */
    sprintf(fn,"%s.lrr.soc",EndoName);
    MESSAGE(1,("Loading socle basis\n"));
    if ((tmp = MatLoad(fn)) == NULL)
	return -1;
    tmp2 = MatInverse(tmp);
    MatFree(tmp);
    tmp = MatTransposed(tmp2);
    MatFree(tmp2);
    head = MatCutRows(tmp,0,headdim);

    return 0;
}




static int Init(int argc, const char **argv)

{
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)                           
        return -1;
    if (ParseArgs() != 0)
        return -1;
    if (ReadFiles() != 0)
        return -1;
    return 0;
}



static void Cleanup()

{
    AppFree(App);
}




static int WriteOutput(Matrix_t *bas)

{
    char name[200];
    MatRep_t *rep = NULL;
    int i;

    /* Write the decomposition basis.
       ------------------------------ */
    sprintf(name,"%s.dec", ModName);
    MESSAGE(1,("Writing the decomposition basis (%s)\n",name));
    MatSave(bas,name);

    /* Transform the generators.
       ------------------------- */
    if (TransformGenerators || WriteAction)
    {
	MESSAGE(1,("Transforming the generators\n"));
	sprintf(name,"%s.std",ModName);
	rep = MrLoad(name,ModInfo.NGen);
	if (rep == NULL)
	{
	    MTX_ERROR("Cannot load generators");
	    return 1;
	}
	if (MrChangeBasis(rep,bas) != 0)
	{
	    MTX_ERROR("Error changing basis");
	    return 1;
	}
	if (TransformGenerators)
	{
	    sprintf(name,"%s.dec",ModName);
	    MESSAGE(1,("Writing transformed generators (%s.1, ...)\n",name));
	    MrSave(rep,name);
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
	    for (k = 0; k < LrrInfo.NCf; k++)
	    {
		int l;
		for (l = 0; l < LrrInfo.Cf[k].dim / LrrInfo.Cf[k].spl; l++)
		{
		    Matrix_t *tmp = MatCut(rep->Gen[i],block_start,block_start,
			compdim[k],compdim[k]);
		    block_start += compdim[k];
		    sprintf(name, "%s.comp%d%c%d.%d", ModName,compdim[k], 
			compnm[k],l+1,i+1);
		    MatSave(tmp, name);
		    MatFree(tmp);
		}
	    }
	}
    }

    return 0;
}


int main(int argc, const char **argv)

{
    int rc = 0;
    int i, j, l, num, dim = 0;
    Matrix_t *bas, *partbas = NULL, *mat, *ker;
    PTR headptr;
    FEL f;
    char name[200];
    FPoly_t *pol = NULL;

                                                                                
    if (Init(argc,argv) != 0)                                                   
    {                                                                           
        MTX_ERROR("Initialization failed");                                     
        return 1;                                                               
    }                                                                           

/* -------------------------------------------------------
   makes the corresponding element of the endomorphismring
   ------------------------------------------------------- */

    bas = MatAlloc(FfOrder,moddim,moddim);
    headptr = head->Data;
    for (i = 0; i < LrrInfo.NCf; i++)
    {
	MESSAGE(1,("Next constituent is %s%s\n",LrrInfo.BaseName,
	    Lat_CfName(&LrrInfo,i)));
	for (j = 0; j < LrrInfo.Cf[i].dim / LrrInfo.Cf[i].spl; j++)
	{
	    num = LrrInfo.Cf[i].dim;
	    do
	    {
		if (partbas != NULL)
		    MatFree(partbas);
		partbas = MatAlloc(FfOrder, moddim, moddim);
		if (num-- == 0)
		{
		    MTX_ERROR("na, most mi van?");
		    return 1;
		}
		for (l = 0; l < enddim; l++)
		{
		    if ((f = FfExtract(headptr, l)) == FF_ZERO)
		    	continue;
		    sprintf(name, "%s.%d", EndoName, l+1);
		    if ((mat = MatLoad(name)) == NULL)
			return 1;
		    MatAddMul(partbas,mat,f);
		    MatFree(mat);
		}
		FfSetNoc(enddim);
		FfStepPtr(&headptr);

		if (pol != NULL)
		    FpFree(pol);
		pol = CharPol(partbas);
	    }
	    while (LrrInfo.Cf[i].dim != 1 && pol->NFactors == 1 
		&& pol->Factor[0]->Degree == 1 && pol->Factor[0]->Data[0] == 0 
		&& pol->Factor[0]->Data[1] == 1); /* i.e.,charpol == x^enddim */
	    FfSetNoc(enddim);
	    headptr = FfGetPtr(headptr,num);


            /* Make the stable kernel.
               ----------------------- */
	    StablePower_(partbas,NULL,&ker);
	    compdim[i] = moddim - ker->Nor;
	    for (l = i - 1; l >= 0 && compdim[i] != compdim[l]; l--)
		;
	    if (l >= 0)
		compnm[i] = (char)(compnm[l] + 1);
	    else
		compnm[i] = 'a';
	    MatFree(ker);
	    MESSAGE(0,("The %d-th direct summand is: %d%c\n\n", j,
		compdim[i], compnm[i]));


	    /* Append <partbas> to <bas>.
	       -------------------------- */
	    MatEchelonize(partbas);
	    MatCopyRegion(bas,dim,0,partbas,0,0,-1,-1);
	    dim += partbas->Nor;
	}
    }
    MatFree(partbas);

    if (dim != moddim)
    {
	MTX_ERROR("es most mi van?");
	for (i = 0; i < LrrInfo.NCf; i++)
    	    printf("%d  ", compdim[i]);
    	printf("\n");
	return 1;
    }



    rc = WriteOutput(bas);


    Cleanup();
    return rc;
}



/**
@page prog_decomp decomp - Decompose a Module                                                  

<<<<<<< HEAD
@section syntax Command Line
=======
@section decomp_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
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

<<<<<<< HEAD
@section inp Input Files
=======
@section decomp_inp Input Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
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

<<<<<<< HEAD
@section out Output Files
=======
@section decomp_out Output Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
@par @em Module.dec
  Basis of the module reflecting the decomposition.
@par @em Module.dec.1, @em Module.dec.2, ...
  Basis of the module reflecting the decomposition.
  Generators in decomp. basis (with -t).
@par @em Module.@em Comp.1, @em Module.@em Comp.2, ...
  Generators on the components (with -a).

<<<<<<< HEAD
@section desc Description
=======
@section decomp_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
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
                                                                                
<<<<<<< HEAD
@section impl Implementation Details
=======
@section decomp_impl Implementation Details
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
  The algorithm used by this program was developed by Magdolna Szöke [@ref Sz98].

**/
