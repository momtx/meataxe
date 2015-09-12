////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Chop with known peak words
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>

#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>



MTX_DEFINE_FILE_INFO
static const MtxApplicationInfo_t AppInfo = {
"pseudochop", "Chop by peakwords",
"SYNTAX\n"
"    pseudochop [-s] <mod> <ref>\n"
"\n"
"ARGUMENTS\n"
"    -s ........... Assume modules are semisimple.\n"
"    <mod> ........ Module to chop.\n"
"    <ref> ........ Reference module (chop and pwkond -t have been run).\n"
"\n"
"FILES\n"
"    <ref>1a.1 ...             i  all compositionfactors occuring in <<gen>.j>\n"
"    <ref>.cfinfo              i  constituent information\n"
"    <mod>.1 ... <mod>.nbgen   i  generators of representation to be chopped\n"
"    <mod>1a.std.1             o  copies of <mod>1a.std.1....\n"
"    <mod>1a.op ...            o  copies of <mod>1a.op ...\n"
"    <mod>1a.k ...             o  see pwkond ...\n"
"    <mod>1a.1 ...             o  copies of <mod>1a.1 ...\n"
"    <mod>.cfinfo              o  constituent information\n"
   };
static MtxApplication_t *App = NULL;
static int semisimp =0;

#define MAX_NAME 500
#define MAX_NBGEN 100

static int ParseCommandLine()
{
    semisimp = AppGetOption(App,"-s --assume-semisimple");
    if (AppGetArguments(App,2,2) != 2)
	return -1;
    if (semisimp)
        printf("Assuming that the representation is semisimple.\n");
    return 0;
}

/*----------------------------------------------------------------------------*/

int main(int argc, const char *argv[])
{
    int dim = 0, i, j;
    char name[MAX_NAME], name2[MAX_NAME];
    MatRep_t *gens;
    Matrix_t *old, *mat;
    Matrix_t *nulsp;
    WgData_t *rep;
    Lat_Info mycfinfo;
    IntMatrix_t *OpTable;
    MtxInitLibrary();
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return 1;
    if (ParseCommandLine() != 0)
	return 1;

    /* Read <ref>.cfinfo */
    if (Lat_ReadInfo(&mycfinfo,App->ArgV[1]) != 0)
	return -1;

    /* Read generators of <mod>, and set up the word generator */
    gens = MrLoad(App->ArgV[0],mycfinfo.NGen);
    rep = WgAlloc(gens);

    /* Run through all possible constituents and calculate their multiplicity in <mod> */
    for ( j = 0; j < mycfinfo.NCf; j++ )
    {
	int oldnul, newnul;
        if (mycfinfo.Cf[j].peakword == 0)
        {
	    MTX_ERROR("0 is definitly not a peakword! - Did you run mkpeak?");
            return 1;
        }
        old = MatInsert(WgMakeWord(rep,mycfinfo.Cf[j].peakword),mycfinfo.Cf[j].peakpol);
	nulsp = MatNullSpace(old);
	newnul = nulsp->Nor; 
        oldnul = 0;

	/* Find stable nullity */
        while (!semisimp && newnul > oldnul)
        {
	    Matrix_t *newmat;
            oldnul = newnul;
            newmat = MatDup(old);
            MatMul(newmat,old);
            MatFree(old);
            MatFree(nulsp);
            old = MatDup(newmat);
            nulsp = MatNullSpace__(newmat); 
            newnul= nulsp->Nor;
        }
        MatFree(old);
        mycfinfo.Cf[j].mult = newnul / mycfinfo.Cf[j].spl;
        dim += mycfinfo.Cf[j].dim * mycfinfo.Cf[j].mult;

        MESSAGE(0,("%s%s occurs %ld times (total dimension now %d out of %d)\n",
	    App->ArgV[1],Lat_CfName(&mycfinfo,j),mycfinfo.Cf[j].mult,
	    dim,gens->Gen[0]->Nor));

	/* Copy generators, std basis, and .op file for this constituent.
           Note: we do this even if this constituent does not occcur in <mod> */
        sprintf(name, "%s%s.k", App->ArgV[0], Lat_CfName(&mycfinfo,j));
        MatSave(nulsp,name);
        sprintf(name, "%s%s.op", App->ArgV[0], Lat_CfName(&mycfinfo,j));
        sprintf(name2, "%s%s.op", App->ArgV[1], Lat_CfName(&mycfinfo,j));
        OpTable = ImatLoad(name2);
        ImatSave(OpTable,name);
        for ( i = 0; i < mycfinfo.NGen; i++ )
        {
            sprintf(name, "%s%s.%d", App->ArgV[0], Lat_CfName(&mycfinfo,j), i+1);
            sprintf(name2, "%s%s.%d", App->ArgV[1], Lat_CfName(&mycfinfo,j), i+1);
            mat = MatLoad(name2);
            if (mat->Field != gens->Gen[0]->Field)
            {
		MTX_ERROR2("%s: %E",name2,MTX_ERR_INCOMPAT); 
                return -1;
            }
            MatSave(mat, name);
            MatFree(mat);
            sprintf(name, "%s%s.std.%d", App->ArgV[0], Lat_CfName(&mycfinfo,j), i+1);
            sprintf(name2, "%s%s.std.%d", App->ArgV[1], Lat_CfName(&mycfinfo,j), i+1);
            mat = MatLoad(name2);
            MatSave(mat, name);
            MatFree(mat);
        } 
    }
    if (dim < gens->Gen[0]->Nor)
        fprintf(stderr, "The given compositionfactors form only %d from whole dimension %d!\n\n",
		dim, gens->Gen[0]->Nor);

    strcpy(mycfinfo.BaseName,App->ArgV[0]);
    Lat_WriteInfo(&mycfinfo);

    return 0;	
}
    

/**
@page prog_pseudochop pseudochop - Chop with known peak words

@see
- @ref prog_chop
- @ref prog_pwkond

@section pseudochop_syntax Syntax
<pre>
pseudochop [-QVs] @em Module @em Reference
</pre>

@par -Q
Quiet, no messages.
@par -V
Verbose, more messages.
@par -s
Assume modules are semisimple
@par @em Module
Name of the module to chop.
@par @em Reference
Name of the reference module, where chop and pwknd -t have been run.

@section pseudochop_inp Input files
@par @em Reference.cfinfo
Constituent information for the reference module.

@par @em ReferenceNx.1 ...
Composition factors of the reference module.

@par @em Module .1 @em Module .2 ...
Action of the generators on @em Module

@section pseudochop_out Output files
The output is the same as for chop.

@section pseudochop_desc Description
This program can be used to produce the chop output for a given module without
actually doing all the work. To do so, there must be another module (the "reference
module") which contains all constituents that occur in the module to be chopped,
and the reference module must have been chopped, and peak words must have been
calculated.

*/
