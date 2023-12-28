////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Chop with known peak words
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>



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
    semisimp = appGetOption(App,"-s --assume-semisimple");
    if (appGetArguments(App,2,2) != 2)
	return -1;
    if (semisimp)
        printf("Assuming that the representation is semisimple.\n");
    return 0;
}

/*----------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    int dim = 0, i, j;
    char name[MAX_NAME], name2[MAX_NAME];
    MatRep_t *gens;
    Matrix_t *old, *mat;
    Matrix_t *nulsp;
    WgData_t *rep;
    IntMatrix_t *OpTable;
    mtxInitLibrary(argv[0]);
    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return 1;
    if (ParseCommandLine() != 0)
	return 1;

    /* Read <ref>.cfinfo */
    LatInfo_t* mycfinfo = latLoad(App->argV[1]);

    /* Read generators of <mod>, and set up the word generator */
    gens = mrLoad(App->argV[0],mycfinfo->NGen);
    rep = wgAlloc(gens);

    /* Run through all possible constituents and calculate their multiplicity in <mod> */
    for ( j = 0; j < mycfinfo->nCf; j++ )
    {
	int oldnul, newnul;
        if (mycfinfo->Cf[j].peakWord == 0)
	    mtxAbort(MTX_HERE,"No peak word defined - Did you run mkpeak?");
        old = matInsert(wgMakeWord(rep,mycfinfo->Cf[j].peakWord),mycfinfo->Cf[j].peakPol);
	nulsp = matNullSpace(old);
	newnul = nulsp->nor; 
        oldnul = 0;

	/* Find stable nullity */
        while (!semisimp && newnul > oldnul)
        {
	    Matrix_t *newmat;
            oldnul = newnul;
            newmat = matDup(old);
            matMul(newmat,old);
            matFree(old);
            matFree(nulsp);
            old = matDup(newmat);
            nulsp = matNullSpace__(newmat); 
            newnul= nulsp->nor;
        }
        matFree(old);
        mycfinfo->Cf[j].mult = newnul / mycfinfo->Cf[j].spl;
        dim += mycfinfo->Cf[j].dim * mycfinfo->Cf[j].mult;

        MTX_LOGI("%s%s occurs %d times (total dimension now %d out of %d)\n",
	    App->argV[1],latCfName(mycfinfo,j),mycfinfo->Cf[j].mult, dim,gens->Gen[0]->nor);

	/* Copy generators, std basis, and .op file for this constituent.
           Note: we do this even if this constituent does not occcur in <mod> */
        sprintf(name, "%s%s.k", App->argV[0], latCfName(mycfinfo,j));
        matSave(nulsp,name);
        sprintf(name, "%s%s.op", App->argV[0], latCfName(mycfinfo,j));
        sprintf(name2, "%s%s.op", App->argV[1], latCfName(mycfinfo,j));
        OpTable = imatLoad(name2);
        imatSave(OpTable,name);
        for ( i = 0; i < mycfinfo->NGen; i++ )
        {
            sprintf(name, "%s%s.%d", App->argV[0], latCfName(mycfinfo,j), i+1);
            sprintf(name2, "%s%s.%d", App->argV[1], latCfName(mycfinfo,j), i+1);
            mat = matLoad(name2);
            if (mat->field != gens->Gen[0]->field)
            {
		mtxAbort(MTX_HERE,"%s: %s",name2,MTX_ERR_INCOMPAT); 
                return -1;
            }
            matSave(mat, name);
            matFree(mat);
            sprintf(name, "%s%s.std.%d", App->argV[0], latCfName(mycfinfo,j), i+1);
            sprintf(name2, "%s%s.std.%d", App->argV[1], latCfName(mycfinfo,j), i+1);
            mat = matLoad(name2);
            matSave(mat, name);
            matFree(mat);
        } 
    }
    if (dim < gens->Gen[0]->nor)
        fprintf(stderr, "The given compositionfactors form only %d from whole dimension %d!\n\n",
		dim, gens->Gen[0]->nor);

    strcpy(mycfinfo->BaseName,App->argV[0]);
    latSave(mycfinfo);
    latDestroy(mycfinfo);

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
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
