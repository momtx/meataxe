////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Order of a matrix or permutations.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>
#include <string.h>

#define MAXORDER 100000         /* Maximal order */
#define MAXORDER_C 1000000         /* Maximal order on cyclic subspaces */

/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static MtxApplicationInfo_t AppInfo = {
   "zor", "Order of a Matrix or Permutation",
   "SYNTAX\n"
   "    zor [-m <MaxOrder>] [-GQVq] <File>\n"
   "\n"
   "OPTIONS\n"
   MTX_COMMON_OPTIONS_DESCRIPTION
   "    -m ...................... Set highest possible order\n"
   "    -q ...................... Quick mode: Find a lower bound for the order\n"
   "    -G ...................... GAP output\n"
   "\n"
   "FILES\n"
   "    <File> .................. I  A matrix or permutation\n"
};

static MtxApplication_t* App = NULL;
static const char* fileName;
static MtxFile_t* file;
static int maxord = -1; /* Limit set with -m option */
static int opt_q = 0;           /* Quick mode */
static int opt_G = 0;           /* GAP output */

////////////////////////////////////////////////////////////////////////////////////////////////////

static void calculateMatrixOrder()
{
   PTR m1, v;
   PTR base, bend;
   uint32_t dim;
   uint32_t* piv;
   char* ispiv;
   int ord;

   ffSetField(file->header[0]);
   const uint32_t nor = file->header[1];
   const uint32_t noc = file->header[2];
   if (nor != noc) {
      mtxAbort(MTX_HERE,"%s: %s",fileName,MTX_ERR_NOTSQUARE);
   }
   m1 = ffAlloc(noc, noc);
   base = ffAlloc(noc + 1, noc);
   piv = NALLOC(uint32_t,noc);
   ispiv = NALLOC(char,noc);
   memset(ispiv,0,noc);
   v = ffAlloc(1, noc);
   ffReadRows(file,m1,noc,noc);
   ord = 1;
   bend = base;

   for (dim = 0; dim < noc;) {
      PTR start = bend;
      int tord = 0;
      int closed = 0;

      /* Find the next seed vector
         ------------------------- */
      uint32_t i = 0;
      while (i < noc && ispiv[i]) ++i;
      MTX_ASSERT(i < noc);
      ffMulRow(bend,FF_ZERO, noc);
      ffInsert(bend,i,FF_ONE);

      /* Calculate order on the cyclic subspace
         -------------------------------------- */
      do {

         /* Save the vector and extend the basis,
            if the vector is linearly independent.
            -------------------------------------- */
         ffCopyRow(v,bend, noc);
         if (!closed) {
            PTR b = base;
            for (i = 0; i < dim; ++i) {
               FEL f = ffExtract(bend,piv[i]);
               if (f != FF_ZERO) {
                  ffAddMulRow(bend,b,ffNeg(ffDiv(f,ffExtract(b,piv[i]))),noc);
               }
               ffStepPtr(&b,noc);
            }
            FEL f;
            uint32_t pv = ffFindPivot(bend,&f, noc);
            if (pv != MTX_NVAL) {
               piv[dim++] = pv;
               ispiv[pv] = 1;
               ffStepPtr(&bend, noc);
            } else {
               closed = 1;
            }
         }

         /* Apply the matrix.
            ----------------- */
         if (++tord > MAXORDER_C) {
            mtxAbort(MTX_HERE,"zor: Partial order is over %d",MAXORDER_C);
         }
         ffMapRow(bend, v,m1,noc,noc);
      }
      while (ffCmpRows(bend,start,noc));

      // Calculate l.c.m. of all partial orders.
      for (i = ord; ord % tord != 0; ord += i) {
      }
      if (ord > MAXORDER) {
         mtxAbort(MTX_HERE,"zor: Order is over %d",MAXORDER);
      }
      if (opt_q && dim > noc / 10) {break;}
      if (maxord > 1) {
         if (ord > maxord) {
            mtxAbort(MTX_HERE, "zor: Order is over %d",maxord);
         }
         if (ord == maxord) {break;}
      }
   }

   if (opt_q && ord != maxord) {
      MTX_LOGI("ORDER IS A MULTIPLE OF %d",ord);
   } else {
      if (opt_G) {
         printf("MeatAxe.Order := %d;\n",ord);
      } else {
         printf("ORDER IS %d\n",ord);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void calculatePermutationOrder()
{
   const uint32_t numberOfPermutations = file->header[2];

   if (opt_G) {
      printf("MeatAxe.Orders := [");
   }
   for (uint32_t i = 1; i <= numberOfPermutations; ++i) {
      Perm_t* perm = permReadData(file);
      uint32_t order = permOrder(perm);
      permFree(perm);
      if (order < 0) {
         mtxAbort(MTX_HERE,"Error calculating order");
      }
      if (opt_G) {
         printf("%s%lu",i > 1 ? "," : "", (unsigned long)order);
      } else {
         printf("ELEMENT %lu HAS ORDER %lu\n",(unsigned long) i, (unsigned long)order);
      }
   }
   if (opt_G) { printf("];\n");}
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char** argv)
{
   // Process command line options.
   App = appAlloc(&AppInfo,argc,argv);
   opt_G = appGetOption(App,"-G --gap");
//   if (opt_G) {
//      MtxMessageLevel = -100;
//   }
   opt_q = appGetOption(App,"-q --quick");
   maxord = appGetIntOption(App,"-m --max-order",-1,1,1000000);

   // Process arguments.
   appGetArguments(App,1,1);
   fileName = App->argV[0];
}


////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
   init(argc, argv);
   file = mfOpen(fileName, "rb");
   mfReadHeader(file);
   uint32_t objectType = mfObjectType(file);
   switch (objectType) {
      case MTX_TYPE_MATRIX:
         calculateMatrixOrder();
         break;
      case MTX_TYPE_PERMUTATION:
         calculatePermutationOrder();
         break;
      default:
         mtxAbort(MTX_HERE,
            "%s: unsupported object type 0x%lx",
            fileName,
            (unsigned long) objectType);
   }
   mfClose(file);
   appFree(App);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zor zor - Order

@section syntax Command Line
<pre>
zor @em Options [-q] [-m @em MaxOrder] @em Input
</pre>

@par @em Options
Standard options, see @ref prog_stdopts
@par -q
Quick mode (see description).
@par -m @em MaxOrder
Set upper limit for the oder of cyclic subspaces (see description).
@par @em Input
Input file.

@section zor_inp Input Files
@par @em Input
Matrix or permutation.

@section zor_desc Description
This program reads a file, containing either permutations, or a
square matrix, and calculates the order(s) and prints the message
<pre>
ORDER IS xxxx</pre>
There are two options to reduce the run time of the program.
Using the "-m" option you can specify a maximal expected order.
If, during the algorithm described below, the order reaches
this limit, the program will stop and print an appropriate message.
The second option, "-q", makes @b zor stop if the dimension of W
(see below) reaches 1/10 of the dimension of the whole space. In
this case, the message is
<pre>
ORDER IS A MULTIPLE OF @em NNN
</pre>
Note: The "-q" and "-m" options have no effect for permutations.

@section zor_impl Implementation Details
If the input is a matrix, the order is found by calculating the orders on cyclic
subspaces and taking the least common multiple.
The algorithm works as follows:
- Let A be the given matrix and V the space A acts upon.
Set W:={0} (the trivial subspace) and o:=1.
- (NEXT) Chose a vector v not in W.  Calculate the cyclic subspace
C generated by v and the order o' of A on C.
- Set o:=lcm(o,o')
- W:=W + C. If W=V, o is the order of A and the program terminates.
Otherwise, continue with (NEXT).

Gaussian elimination is used to maintain a basis of W in echelon form.
In order to avoid infinite loops, there is a limit on o'. If the
vector does not return after 1000 multiplications the order is assumed
to be infinite and the program stops with an error message. This happens
also if the value of o exceeds 100000.

If the input file contains permutations, each one is read in and
its order is calculated as the least common multiple of the orbit
sizes. The result is printed in the format
<pre>
ELEMENT nn HAS ORDER nnn
</pre>

The whole matrix plus a second matrix of the same size must fit into
memory. In the case of permutations, there must be enough memory for
one permutation.
 **/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
