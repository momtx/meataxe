/* ============================= C MeatAxe ==================================
   File:        $Id: zmw.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Make word.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"

#include <ctype.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

static int opt_G = 0;
static int NGen = -1;
static Matrix_t *Gen[MAXGEN];
static int WordNo, WordNo2 = -1;
static MatRep_t *Rep = NULL;
static WgData_t *WGen = NULL;
static Poly_t *Poly = NULL;
static const char *WordFileName = NULL;
static const char *NspFileName = NULL;


static MtxApplicationInfo_t AppInfo = { 
"zmw", "Make word", 
"SYNTAX\n"
"    zmw [<Options>] <No> <Gen1> <Gen2> [<Word> [<Nsp>]]\n"
"    zmw [<Options>] -g <NGen> <No> <Gen> [<Word> [<Nsp>]]\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... GAP output\n"
"    -g ...................... Set number of generators\n"
"\n"
"ARGUMENTS\n"
"    <No> .................... Word or word range (A-B), optionally followed\n"
"                              by polynomial (A/c,c,c...c)\n"
"    <Gen1>,<Gen2> ........... Generators\n"
"    <Gen> ................... Generator base name (with -g)\n"
"    <Word> .................. Output file name\n"
"    <Nsp> ................... Null-Space file name\n"
"\n"
"FILES\n"
"    <Gen>.{1,2...} .......... I Generators (with -g)\n"
"    <Gen1>, <Gen2> .......... I Generators (without -g)\n"
"    <Word> .................. O Word\n"
"    <Nsp> ................... O Null-Space of <Word>\n"

};

static MtxApplication_t *App = NULL;





static int ReadGenerators()

{
    int ngen = NGen == -1 ? 2 : NGen;
    int i;
    char fn0[200];

    for (i = 0; i < ngen; ++i)
    {
	char fn[200];
	if (NGen != -1)
	    sprintf(fn,"%s.%d",App->ArgV[1],i + 1);
	else
	    strcpy(fn,App->ArgV[1 + i]);
	if ((Gen[i] = MatLoad(fn)) == NULL)
	    return -1;
	if (i > 0)
	{
	    if (Gen[i]->Field != Gen[0]->Field 
		|| Gen[i]->Nor != Gen[0]->Nor || Gen[i]->Noc != Gen[0]->Noc)
	    {
		MTX_ERROR3("%s and %s: %E",App->ArgV[1],App->ArgV[i+1],
		    MTX_ERR_INCOMPAT);
	    }
	}
	else
	    strcpy(fn0,fn);
    }

    NGen = ngen;
    if ((Rep = MrAlloc(NGen,Gen,0)) == NULL)
	return -1;
    if ((WGen = WgAlloc(Rep)) == NULL)
	return -1;

    return 0;
}





/* ------------------------------------------------------------------
   ParseWord() - Parse word/polynomial specification

   Examples:
        spec          WordNo   WordNo2   Poly
	------------  ------   -------   ------------
        10            10       -1        NULL
	1-100         1        100       NULL
	30/1,-1       30       -1        x-1
	30-32/1,1,1   30       31        x^2+x+1
   ------------------------------------------------------------------ */

static int ParseWord(const char *spec)

{
    if (!isdigit(*spec))
	return -1;
    WordNo = atoi(spec);
    while (isdigit(*spec)) ++spec;
    if (*spec == '-')
    {
	++spec;
	if (!isdigit(*spec))
	    return -1;
	WordNo2 = atoi(spec);
	while (isdigit(*spec)) 
	    ++spec;
    }
    if (*spec == '/')
    {
	const char *c;
	int deg;
	++spec;

	for (c = spec, deg = 0; *c != 0; ++c)
	{
	    if (*c == ',') 
		++deg;
	}
	Poly = PolAlloc(FfOrder,deg);
	while (deg >= 0)
	{
	    if (!isdigit(*spec) && (*spec != '-' || !isdigit(spec[1])))
		return -1;
	    Poly->Data[deg--] = FfFromInt(atoi(spec));
	    while (*spec == '-' || isdigit(*spec)) 
		++spec;
	    if (*spec == ',')
		++spec;
	}
    
    }
    if (*spec != 0)
	return -1;

    MESSAGE(1,("Word %d..%d, Poly=",WordNo,WordNo2));
    if (MSG1)
    {
	if (Poly != NULL)
	{
	    PolPrint(NULL,Poly);
	    printf("\n");
	}
	else
	    printf("x\n");
    }
    return 0;
}




/* ------------------------------------------------------------------
   init() - Process command line options and arguments
   ------------------------------------------------------------------ */

static int Init(int argc, const char **argv)

{
    int min_num_args;

    /* Process command line options.
       ----------------------------- */
    App = AppAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;
    opt_G = AppGetOption(App,"-G --gap");
    if (opt_G) MtxMessageLevel = -100;
    NGen = AppGetIntOption(App,"-g",-1,1,MAXGEN);

    /* Process arguments.
       ------------------ */
    min_num_args = (NGen == -1) ? 3 : 2;
    if (AppGetArguments(App,min_num_args,min_num_args + 2) < 0)
	return -1;
    if (App->ArgC > min_num_args)
	WordFileName = App->ArgV[min_num_args];
    if (App->ArgC > min_num_args + 1)
	NspFileName = App->ArgV[min_num_args + 1];

    /* Other initialization.
       --------------------- */
    if (ReadGenerators() != 0)
	return -1;
    if (ParseWord(App->ArgV[0]) != 0)
    {
	MTX_ERROR("Invalid word/polynomial specification");
	return -1;
    }

    return 0;
}




static void Cleanup()

{
    WgFree(WGen);
    MrFree(Rep);
    AppFree(App);
}



static int MakeWord()

{
    Matrix_t *w;
    int lastword = WordNo2 >= 0 ? WordNo2 : WordNo;


    if (Poly != NULL && MSG0)
    {
	    PolPrint("Using polynomial p(x)",Poly);
    }
    MESSAGE(0,("Number Nullity Word\n"));
    for (; WordNo <= lastword; ++WordNo)
    {
	w = WgMakeWord(WGen,WordNo);
	if (w == NULL)
	    return -1;
	if (Poly != NULL)
	    MatInsert_(w,Poly);
	MESSAGE(0,("%6d",WordNo));
	if (WordFileName != NULL && WordNo2 == -1)
	    MatSave(w,WordFileName);
	if (WordNo2 != -1 || NspFileName != NULL)
	{
	    Matrix_t *nsp = MatNullSpace_(w,0);
	    if (WordNo2 == -1)
		MatSave(nsp,NspFileName);
	    MESSAGE(0,("%8d",nsp->Nor));
	    MatFree(nsp);
	}
	else
	    MESSAGE(0,("        "));
	MESSAGE(0,(" %s\n",WgSymbolicName(WGen,WordNo)));
	MatFree(w);
    }
    return 0;
}






/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, const char **argv)

{
    int rc;

    if (Init(argc,argv) != 0)
	return 1;
    rc = MakeWord();
    Cleanup();
    return rc;
}



/**
@page prog_zmw zmw - Make Word

<<<<<<< HEAD
@section syntax Command Line
=======
@section zmw_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
<pre>
zmw @em Options @em No @em Gen1 @em Gen2 [@em Word [@em Nsp]]
zmw @em Options -g @em NGen @em No @em Gen [@em Word [@em Nsp]]
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -g @em NGen
  Set the number of generators.
@par @em No
  Word number/polynomial (see description).
@par @em Gen1, @em Gen2
  Generators (without -g).
@par @em Gen
  Generator base name (with -g).
@par @em Word
  Word.
@par @em Nsp
  Null-space.


<<<<<<< HEAD
@section inp Input Files
=======
@section zmw_inp Input Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
@par @em Gen1, @em Gen2
  Generators (without -g).
@par @em Gen.1, @em Gen.2, ...
  Generators (with -g).

<<<<<<< HEAD
@section out Output Files
=======
@section zmw_out Output Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
@par @em Word
  Word.
@par @em Nsp
  Null-space.

<<<<<<< HEAD
@section desc Description
=======
@section zmw_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
@b zmw calculates an element of the algebra generated by a set of matrices.
The word must be specified by its "word number", see @ref wgen.

The first form of the @b zmw command assumes that there are two generators,
and the generators must be given as two files names. The second form uses
the "-g" option to specify the number of generators. In this case, a base
name is given on the command line, and the generator file names are constructed
by appending ".1", ".2", ... to the base name.

The word is written to @em Word, if present. If @em Nsp is present,
the kernel of the word is calculated and written to this file.

@em No may be a single number, or a range of numbers in the form "a-b".
In this case the nullities of all words are calculated and printed,
but no output file is written, even if @em Word and @em Nsp are present.

@b zmw can insert the words into a polynomial before calculating the
null-space. To use this feature, the polynomial must be appended to 
the word number, separated by a "/".
After the "/", all coefficients, including zeroes,
of the polynomial must be specified as a comma-separated list.
For example,
<pre>
zmw 103/1,1,0,-1 gen1 gen2</pre>
calculates word number 103, and inserts the result into the 
polynomial x<sup>3</sup>+x<sup>2</sup>-1.
**/

