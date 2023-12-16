////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Make word.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <ctype.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>

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

////////////////////////////////////////////////////////////////////////////////////////////////////

static void readGenerators()
{
   int ngen = NGen == -1 ? 2 : NGen;
   int i;
   char fn0[200];

   for (i = 0; i < ngen; ++i) {
      char fn[200];
      if (NGen != -1) {
         sprintf(fn, "%s.%d", App->argV[1], i + 1);
      }
      else {
         strcpy(fn, App->argV[1 + i]);
      }
      Gen[i] = matLoad(fn);
      if (i > 0) {
         if (Gen[i]->field != Gen[0]->field
             || Gen[i]->nor != Gen[0]->nor || Gen[i]->noc != Gen[0]->noc) {
            mtxAbort(MTX_HERE, "%s and %s: %s", App->argV[1], App->argV[i + 1],
               MTX_ERR_INCOMPAT);
         }
      }
      else {
         strcpy(fn0, fn);
      }
   }

   NGen = ngen;
   Rep = mrAlloc(NGen, Gen, 0);
   WGen = wgAlloc(Rep);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

//   Parse word/polynomial specification   Examples:
//
//   spec          WordNo   WordNo2   Poly
//   ------------  ------   -------   ------------
//   10            10       -1        NULL
//   1-100         1        100       NULL
//   30/1,-1       30       -1        x-1
//   30-32/1,1,1   30       31        x^2+x+1

static int parseWord(const char* spec)
{
   if (!isdigit(*spec)) {
      return -1;
   }
   WordNo = atoi(spec);
   while (isdigit(*spec)) {
      ++spec;
   }
   if (*spec == '-') {
      ++spec;
      if (!isdigit(*spec)) {
         return -1;
      }
      WordNo2 = atoi(spec);
      while (isdigit(*spec)) {
         ++spec;
      }
   }
   if (*spec == '/') {
      const char* c;
      int deg;
      ++spec;

      for (c = spec, deg = 0; *c != 0; ++c) {
         if (*c == ',') {
            ++deg;
         }
      }
      Poly = polAlloc(ffOrder, deg);
      while (deg >= 0) {
         if (!isdigit(*spec) && (*spec != '-' || !isdigit(spec[1]))) {
            return -1;
         }
         Poly->data[deg--] = ffFromInt(atoi(spec));
         while (*spec == '-' || isdigit(*spec)) {
            ++spec;
         }
         if (*spec == ',') {
            ++spec;
         }
      }

   }
   if (*spec != 0) {
      return -1;
   }

   MTX_XLOGD(msg) {
      sbPrintf(msg, "Word %d..%d, Poly=", WordNo, WordNo2);
      if (Poly != NULL) {
         polFormat(msg, Poly);
      }
      else {
         sbAppend(msg, "x");
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Process command line options and arguments

static int Init(int argc, char** argv)
{
   int min_num_args;

   App = appAlloc(&AppInfo, argc, argv);
   opt_G = appGetOption(App, "-G --gap");
//   if (opt_G) { MtxMessageLevel = -100; }
   NGen = appGetIntOption(App, "-g", -1, 1, MAXGEN);

   min_num_args = (NGen == -1) ? 3 : 2;
   appGetArguments(App, min_num_args, min_num_args + 2);
   if (App->argC > min_num_args) {
      WordFileName = App->argV[min_num_args];
   }
   if (App->argC > min_num_args + 1) {
      NspFileName = App->argV[min_num_args + 1];
   }

   readGenerators();
   if (parseWord(App->argV[0]) != 0) {
      mtxAbort(MTX_HERE, "Invalid word/polynomial specification");
      return -1;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanup()
{
   if (Poly != NULL)
      polFree(Poly);
    wgFree(WGen);
    mrFree(Rep);
    appFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int MakeWord()
{
   Matrix_t* w;
   int lastword = WordNo2 >= 0 ? WordNo2 : WordNo;

   if (Poly != NULL) {
      MTX_XLOGI(sb) {
         sbAppend(sb, "Using polynomial p(x)=");
         polFormat(sb, Poly);
      }
   }
   MTX_LOGI("Number Nullity Word");
   for (; WordNo <= lastword; ++WordNo) {
      w = wgMakeWord(WGen, WordNo);
      if (w == NULL) {
         return -1;
      }
      if (Poly != NULL) {
         matInsert_(w, Poly);
      }
      if (WordFileName != NULL && WordNo2 == -1) {
         matSave(w, WordFileName);
      }
      Matrix_t* nsp = NULL;
      if (WordNo2 != -1 || NspFileName != NULL) {
         nsp = matNullSpace_(w, 0);
         if (WordNo2 == -1) {
            matSave(nsp, NspFileName);
         }
      }

      MTX_XLOGI(msg) {
         sbPrintf(msg, "%6d", WordNo);
         if (nsp != NULL) {
            sbPrintf(msg, "%8d", nsp->nor);
         }
         else {
            sbAppend(msg, "        ");
         }
         sbPrintf(msg, " %s", wgSymbolicName(WGen, WordNo));
      }

      if (nsp) matFree(nsp);
      matFree(w);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    if (Init(argc,argv) != 0)
	return 1;
    int rc = MakeWord();
    cleanup();
    return rc;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zmw zmw - Make Word

@section zmw_syntax Command Line
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


@section zmw_inp Input Files
@par @em Gen1, @em Gen2
  Generators (without -g).
@par @em Gen.1, @em Gen.2, ...
  Generators (with -g).

@section zmw_out Output Files
@par @em Word
  Word.
@par @em Nsp
  Null-space.

@section zmw_desc Description
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

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
