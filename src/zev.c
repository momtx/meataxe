////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculate eigenvalues and multiplitcities.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <ctype.h>
#include <string.h>
#include <stdlib.h>



#define MAXDEG 200		/* Max. degree */
#define LINEWIDTH 2048		/* Line width for input file */


static FILE *src; 		        /* Input file */
static Matrix_t *Matrix;		/* Matrix A */
static const char *matname;
static char grpname[LINEWIDTH] = "";	/* Input selector (empty = ALL) */
static uint32_t deg;			/* Degree */
static char name[LINEWIDTH];		/* Value (atlas notation) */
static char thisgrp[LINEWIDTH];
static int opt_G = 0;			/* -G option (GAP output) */


static MtxApplicationInfo_t AppInfo = { 
"zev", "Eigenvalues",
"SYNTAX\n"
"    zev [-GQV] <Matrix> [<Poly> [<Group>]]\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G, --gap ............... GAP output (implies -Q)\n"
"\n"
"ARGUMENTS\n"
"    <Matrix> ................ A square matrix\n"
"    <Poly> .................. Data file with polynomials (default: standard input)\n"
"    <Group> ................. Selects a group of polynomials (default: all)\n"
"\n"
"FILES\n"
"    <Matrix> ................ I The matrix\n"
"    <Poly> .................. I Polynomial definition file\n"
};



static MtxApplication_t *App = NULL;


////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    opt_G = appGetOption(App,"-G --gap");
//    if (opt_G) MtxMessageLevel = -100;
    switch (appGetArguments(App,1,3)) {
	case 3:
	    strcpy(grpname,App->argV[2]);
	    /* NO BREAK! */
	case 2:
	    if (strcmp(App->argV[1],"-"))
	    {
		src = sysFopen(App->argV[1],"rb");
	    }
	    else
		src = stdin;
	    break;
        case 1:
            src = stdin;
            break;
    }
    matname = App->argV[0];
    Matrix = matLoad(matname);
    if (Matrix->nor != Matrix->noc)
	mtxAbort(MTX_HERE,"%s: %s",matname,MTX_ERR_NOTSQUARE);
    ffSetField(Matrix->field);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Calculate nullity and print the result. Releases the matrix.

void Gauss(Matrix_t* w)
{
    static int first = 1;
    int nullity = matNullity__(w);
    if (opt_G)	/* GAP output */
    {	
	int mult = nullity / deg;

    	if (mult > 0)
    	{
	    if (!first)
		printf(" + ");
	    else
	    {
		printf("MeatAxe.BrauerChar := ");
		first = 0;
	    }
	    printf("%d*(%s)",mult,name);
	}
	if (nullity % deg != 0)
	    fprintf(stderr,"Non-integer multiplicity for %s\n",name);
    }
    else	/* Table output */
    {
	printf("%10s %20s %4ld %4ld",thisgrp,name,(unsigned long)deg,(unsigned long)nullity/deg);
	if (nullity % deg != 0)
	    printf(" (non-integer)\n");
	else
	    printf("\n");
	fflush(stdout);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads one input line. Skips comments (lines starting with '#').
/// Returns 1 on success or 0 at end-of-file.

int readLine(char *buf)
{
   while (fgets(buf,LINEWIDTH,src) != NULL) {
      if (*buf == '#') continue;
      char *c = buf + strlen(buf);
      while (c > buf && (isspace(c[-1]))) --c;
      *c = 0;
      return 1;
   }
   return 0; // end of file
}


////////////////////////////////////////////////////////////////////////////////////////////////////

// Read the next polynomial

Poly_t* getNextPolynomial(void)
{
    char line[LINEWIDTH];
    while (!feof(src))
    {	
	if (!readLine(line)) 
	    return NULL;
	if (*line != ' ')	/* new group */
	{	
	    strcpy(thisgrp,line);
	    continue;
	}
	if (grpname[0] == 0 || !strcmp(thisgrp,grpname))
	    break;
    }
    char* c = strtok(line," \t");
    strcpy(name,c);
    FEL coeff[MAXDEG+1];
    for (deg = 0; (c = strtok(NULL," \t")) != NULL; ++deg)
    {	
        coeff[deg] = (strcmp(c,"-1") == 0) ? ffNeg(FF_ONE) : ffFromInt(atoi(c));
    }
    --deg;
    Poly_t* poly = polAlloc(ffOrder,deg);
    for (uint32_t i = 0; i <= deg; ++i)
	poly->data[i] = coeff[deg - i];
    return poly;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    init(argc,argv);
    Poly_t* poly;
    while ((poly = getNextPolynomial()) != NULL)
    {
	Matrix_t* w = matInsert(Matrix,poly);
        polFree(poly);
	Gauss(w);
    }
    if (opt_G) printf(";\n");
    matFree(Matrix);
    appFree(App);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zev zev - Eigenvalues

@section zev_syntax Command Line
<pre>
zev @em Options [-G] @em Matrix [@em Poly [@em Group]]
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts.
@par -G
  GAP output.
@par @em Matrix
  Input matrix.
@par @em Poly
  Polynomials file name (default: standard input).
@par @em Group
  Group of polynomials to check (default: all groups).

@section zev_inp Input Files
@par @em Matrix
  Input matrix.
@par @em Poly
  Polynomials.


@section zev_desc Description
This program reads a matrix from @em Matrix and a list of polynomials
from @em Poly (or from the standard input).
For each input polynomial, it inserts the matrix into the polynomial,
calculates the nullity, and puts out this nullity, divided by the degree,
along with a text from the input.

The program was specifically designed to assist in the calculation of
the Brauer characters of diagonalizable matrices, with the text giving
the complex number which is the Brauer character of the companion
matrix for that polynomial. Usually the polynomials have been prepared
in a separate data file and are fed into @b zev by giving the file name or
by redirecting its input. The preparation of the input polynomials
is generally a time-consuming task if it is done by hand, but there are
data files available for the most commonly used fields. These
files should be located in the library directory.
They are distributed with this release of the C
MeatAxe. If the user is familiar with the computer program system
GAP, he will find it easy to create his own data files.

If the nullity is not a multiple of the degree, @b zev prints a
warning message.


@subsection pff Polynomial File Format
The data file contains the polynomials in text form.
Several polynomials can be comprised in a group, and the data file
can contain any number of groups of polynomials.
This allows several sets of polynomials to be kept in one data file
(for example, all polynomials for a given field), the appropriate
polynomials being selected trough the @em Group argument on the
command line.

The file is read and interpreted line by line. There are three
types of lines:
- Comment lines, beginning with a "#". These lines
  are simply ignored by @b zev, as are empty lines.
- Group headers. Each line beginning with a non-space character
  is interpreted as the beginning of a new group of polynomials.
  Such lines contain only one text field, the name of the group
  (up to 1023 characters).
- Lines beginning with a space are interpreted as polynomials.
  The format is
<pre>
[space]@em Name @em Coefficients
</pre>
  where @em Name is any text (up to 1023 characters), and
  @em Coefficients are the coefficients of the polynomial (in
  free format). Note that the first character must be an ordinary
  space charcter, a TAB is not allowed!
  The coefficients must use the names as specified by the arithmetic
  --- usually 0,1,...q-1.  The one exception is -1, which
  the program treats specially as "0-1" so that the cyclotomic
  polynomials can be used over all fields. The coefficients
  are in decreasing degree, starting with the coefficient of the
  highest power of x and continuing, ending up with the constant
  term.

Here is an example:
<pre>
# Sample input file for zev
# Some polynomials over GF(5)
#
p11b11
  1         1 4
  b11       1 4 4 1 3 4
  -1-b11    1 2 4 1 1 4
p13c13
  1         1 4
  c13       1 3 0 3 1
  c13*3     1 1 4 1 1
  c13*9     1 2 1 2 1
</pre>
This file contains 7 polynomials in two groups. The polynomial
"b11" in group "p11b11" is x<sup>5</sup>+4x<sup>4</sup>+4x<sup>3</sup>+x<sup>2</sup>+3x+4.


@subsection zev_output_format Output Format
There are two output formats. By default the nullities are printed in
tabular form giving group, name, degree and multiplicity (i.e., nullity
divided by degree) for each
polynomial. If the "-G" option is given, @b zev prints an algebraic
expression which can be read from GAP.
Here is an example with an 8 by 8 matrix over GF(3), polynomials being
read from "poly3":
<pre>
$ zev mat poly3 p8i2
     p8i2                    1    1    1
     p8i2                   -1    1    4
     p8i2                    0    2    2
     p8i2                   i2    2    0
     p8i2                  -i2    2    1
$ zev -G mat poly3 p8i2
MeatAxe.BrauerChar := 1*(1) + 4*(-1) + 2*(0) + 1*(-i2);
</pre>
Note that "i2" does not appear in the expression because its
coefficient is zero.

@section zev_impl Implementation Details
There must be enough memory to hold the input matrix and two
more matrices of the same size.
Lines in the polynomial input file must not be longer than 1023
characters. 

It is not checked that the input file is a matrix.
TAB characters at the beginning of a line are not interpreted
as space.

*/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
