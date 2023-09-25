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

static MtxApplication_t* App = NULL;

static int opt_G = 0;           // GAP output */
static const char* fileNameInp;
static const char* fileNameOut;
static MtxFile_t* fileOut = NULL;
static enum {M_E2, M_E3, M_E4, M_S2} mode;
static uint32_t objectType = 0;
static uint32_t field = 0;              // Field
static uint32_t nor, noc;               // Input size
static uint32_t norOut, nocOut;         // Output size

static Matrix_t* matrixInp = NULL;
static PTR* rowIn = NULL;               // Pointers to the rows of xInt
static PTR rowOut = NULL;               // One row of the output matrix
static Perm_t* permInp;
static Perm_t* permOut;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void prepare()
{
   // Read input file.
   MtxFile_t* f = mfOpen(fileNameInp);
   mfReadHeader(f);
   objectType = mfObjectType(f);
   if (mode != M_E2 && mode != M_S2 && mode != M_E3 && objectType != MTX_TYPE_MATRIX) {
      mtxAbort(MTX_HERE, "%s: %s", fileNameInp, MTX_ERR_NOTMATRIX);
   }
   if (objectType == MTX_TYPE_MATRIX) {
      field = f->header[0];
      nor = f->header[1];
      noc = f->header[2];
      ffSetField(field);
      matrixInp = matReadData(f->file, f->header);
      rowIn = NALLOC(PTR, nor);
      for (uint32_t i = 0; i < nor; ++i) {
         rowIn[i] = matGetPtr(matrixInp, i);
      }
   } else if (objectType == MTX_TYPE_PERMUTATION) {
      permInp = permReadData(f->file, f->header);
      field = 0;
      nor = noc = permInp->degree;
   } else {
      mtxAbort(MTX_HERE,
               "%s: unsupported object type 0x%lx",
               fileNameInp,
               (unsigned long) objectType);
   }
   mfClose(f);

   // Calculate output size.
   uint64_t nor2_ = 0;
   uint64_t noc2_ = 0;
   switch (mode) {
      case M_S2:
         nor2_ = ((uint64_t)nor * (nor + 1)) / 2;
         noc2_ = field >= 2 ? ((uint64_t)noc * (noc + 1)) / 2 : nor2_;
         break;
      case M_E2:
         MTX_ASSERT(nor > 1);
         nor2_ = ((uint64_t)nor * (nor - 1)) / 2;
         if (field >= 2) {
            MTX_ASSERT(noc > 0);
            noc2_ = (uint64_t)noc * (noc - 1) / 2;
         } else {
            noc2_ = nor2_;
         }
         break;
      case M_E3:
         MTX_ASSERT(nor > 2);
         nor2_ = (uint64_t)nor * (nor - 1) / 2 * (nor - 2) / 3;
         if (field >= 2) {
            MTX_ASSERT(noc > 2);
            noc2_ = (uint64_t)noc * (noc - 1) / 2 * (noc - 2) / 3;
         } else {
            noc2_ = nor2_;
         }
         break;
      case M_E4:
         MTX_ASSERT(nor > 3);
         nor2_ = (uint64_t)nor * (nor - 1) / 2 * (nor - 2) / 3 * (nor - 3) / 4;
         if (objectType == MTX_TYPE_MATRIX) {
            MTX_ASSERT(noc > 3);
            noc2_ = (uint64_t)noc * (noc - 1) / 2 * (noc - 2) / 3 * (noc - 3) / 4;
         } else {
            noc2_ = nor2_;
         }
         break;
      default:
         mtxAbort(MTX_HERE, "Unknown mode %d", (int)mode);
   }
   // Check for overflow
   MTX_ASSERT((nor2_ >> 32) == 0);
   MTX_ASSERT((noc2_ >> 32) == 0);
   norOut = (uint32_t) nor2_;
   nocOut = (uint32_t) noc2_;

   // Prepare output buffer and file.
   if (objectType == MTX_TYPE_MATRIX) {
      MESSAGE(0, ("Output is %ld x %ld\n", (long)norOut, (long)nocOut));
      rowOut = ffAlloc(1, nocOut);
      fileOut = mfCreate(fileNameOut, field, norOut, (field >= 2) ? nocOut : 1);
   } else {
      MESSAGE(0, ("Output has degree %ld\n", (long)norOut));
      fflush(stdout);
      permOut = permAlloc(norOut);
   }
}


/* ------------------------------------------------------------------
   zs2() - Symmetric square
   ------------------------------------------------------------------ */

static void zs2()
{
   long i1, i2, j1, j2, j3;
   FEL f11, f12, f21, f22;
   FEL w1, w2, f1, f2;

   MESSAGE(1, ("Mode S2, part 1\n"));
   for (i1 = 0; i1 < nor - 1; ++i1) {
      for (i2 = i1 + 1; i2 < nor; ++i2) {
         ffMulRow(rowOut, FF_ZERO, nocOut);
         j3 = 0;
         for (j1 = 0; j1 < noc - 1; ++j1) {
            f11 = ffExtract(rowIn[i1], j1);
            f21 = ffExtract(rowIn[i2], j1);
            for (j2 = j1 + 1; j2 < noc; ++j2) {
               f12 = ffExtract(rowIn[i1], j2);
               f22 = ffExtract(rowIn[i2], j2);
               w1 = ffMul(f11, f22);
               w2 = ffMul(f12, f21);
               ffInsert(rowOut, j3, ffAdd(w1, w2));
               ++j3;
            }
         }
         for (j2 = 0; j2 < noc; ++j2) {
            f1 = ffExtract(rowIn[i1], j2);
            f2 = ffExtract(rowIn[i2], j2);
            ffInsert(rowOut, j3, ffMul(f1, f2));
            ++j3;
         }
         mfWriteRows(fileOut, rowOut, 1, nocOut);
      }
   }

   MESSAGE(1, ("Mode S2, part 2\n"));
   for (i1 = 0; i1 < nor; ++i1) {
      j3 = 0;
      ffMulRow(rowOut, FF_ZERO, nocOut);
      for (j1 = 0; j1 < noc - 1; ++j1) {
         f1 = ffExtract(rowIn[i1], j1);
         for (j2 = j1 + 1; j2 < noc; ++j2) {
            f2 = ffExtract(rowIn[i1], j2);
            w2 = ffMul(f1, f2);
            ffInsert(rowOut, j3, ffAdd(w2, w2));
            ++j3;
         }
      }
      for (j2 = 0; j2 < noc; ++j2) {
         f1 = ffExtract(rowIn[i1], j2);
         ffInsert(rowOut, j3, ffMul(f1, f1));
         ++j3;
      }
      mfWriteRows(fileOut, rowOut, 1, nocOut);
   }
}


/* ------------------------------------------------------------------
   zs2p() - Antisymmetric square (permutations)
   ------------------------------------------------------------------ */

static uint32_t maps2(uint32_t i, uint32_t k)
{
   if (i <= k) {
      return (k * (k + 1)) / 2 + i;
   } else {
      return (i * (i + 1)) / 2 + k;
   }
}


static void zs2p()
{
   const uint32_t* p1 = permInp->data;
   uint32_t* p2 = permOut->data;
   int i;

   for (i = 0; i < nor; ++i) {
      int k;
      for (k = 0; k <= i; ++k) {
         p2[maps2(i, k)] = maps2(p1[i], p1[k]);
      }
   }
   permSave(permOut, fileNameOut);
}


/* ------------------------------------------------------------------
   ze2p() - Antisymmetric square (permutations)
   ------------------------------------------------------------------ */

static uint32_t mape2(uint32_t i, uint32_t k)
{
   if (i < k) {
      return (k * (k - 1)) / 2 + i;
   } else {
      return (i * (i - 1)) / 2 + k;
   }
}


static void ze2p()
{
   const uint32_t* p1 = permInp->data;
   uint32_t* p2 = permOut->data;

   for (uint32_t i = 0; i < nor; ++i) {
      for (uint32_t k = 0; k < i; ++k) {
         p2[mape2(i, k)] = mape2(p1[i], p1[k]);
      }
   }
   permSave(permOut, fileNameOut);
}


/* ------------------------------------------------------------------
   ze2() - Antisymmetric square
   ------------------------------------------------------------------ */

static void ze2()
{
   for (uint32_t i1 = 0; i1 < nor - 1; ++i1) {
      for (uint32_t i2 = i1 + 1; i2 < nor; ++i2) {
         ffMulRow(rowOut, FF_ZERO, nocOut);
         uint32_t j3 = 0;
         for (uint32_t j1 = 0; j1 < noc - 1; ++j1) {
            const FEL f11 = ffExtract(rowIn[i1], j1);
            const FEL f21 = ffExtract(rowIn[i2], j1);
            for (uint32_t j2 = j1 + 1; j2 < noc; ++j2) {
               const FEL f12 = ffExtract(rowIn[i1], j2);
               const FEL f22 = ffExtract(rowIn[i2], j2);
               const FEL w1 = ffMul(f11, f22);
               const FEL w2 = ffMul(f12, f21);
               const FEL w3 = ffSub(w1, w2);
               ffInsert(rowOut, j3, w3);
               ++j3;
            }
         }
         mfWriteRows(fileOut, rowOut, 1, nocOut);
      }
   }
}


/* ------------------------------------------------------------------
   ze3() - Antisymmetric cube
   ------------------------------------------------------------------ */

static void ze3()
{
   FEL f11, f12, f13, f21, f22, f23, f31, f32, f33;
   FEL e, g12, g13, g23;
   int i1, i2, i3, j1, j2, j3, jins;

   for (i1 = 0; i1 < nor - 2; ++i1) {
      MESSAGE(1, ("i1 = %d\n", i1));
      for (i2 = i1 + 1; i2 < nor - 1; ++i2) {
         MESSAGE(2, ("i2 = %d\n", i2));
         for (i3 = i2 + 1; i3 < nor; ++i3) {
            MESSAGE(3, ("i3 = %d\n", i3));
            ffMulRow(rowOut, FF_ZERO, nocOut);
            jins = 0;
            for (j1 = 0; j1 < noc - 2; ++j1) {
               f11 = ffExtract(rowIn[i1], j1);
               f21 = ffExtract(rowIn[i2], j1);
               f31 = ffExtract(rowIn[i3], j1);
               for (j2 = j1 + 1; j2 < noc - 1; ++j2) {
                  f12 = ffExtract(rowIn[i1], j2);
                  f22 = ffExtract(rowIn[i2], j2);
                  f32 = ffExtract(rowIn[i3], j2);
                  g12 = ffSub(ffMul(f11, f22), ffMul(f21, f12));
                  g13 = ffSub(ffMul(f31, f12), ffMul(f11, f32));
                  g23 = ffSub(ffMul(f21, f32), ffMul(f31, f22));
                  for (j3 = j2 + 1; j3 < noc; ++j3) {
                     f13 = ffExtract(rowIn[i1], j3);
                     f23 = ffExtract(rowIn[i2], j3);
                     f33 = ffExtract(rowIn[i3], j3);
                     e = ffAdd(ffAdd(ffMul(g12, f33), ffMul(g13, f23)),
                               ffMul(g23, f13));
                     ffInsert(rowOut, jins, e);
                     ++jins;
                  }
               }
            }
            mfWriteRows(fileOut, rowOut, 1l, nocOut);
         }
      }
   }
}


/* ------------------------------------------------------------------
   ze3p() - Antisymmetric cube (permutations)
   ------------------------------------------------------------------ */

#define SWAP(x, y) {tmp = x; x = y; y = tmp;}

static uint32_t mape3(uint32_t i, uint32_t k, uint32_t l)
{
   register long tmp;
   if (i < k) { SWAP(i, k);}
   if (i < l) { SWAP(i, l);}
   if (k < l) { SWAP(k, l);}
   return i * (i - 1) / 2 * (i - 2) / 3 + k * (k - 1) / 2 + l;
}


static void ze3p()
{
   uint32_t* p1 = permInp->data;
   uint32_t* p2 = permOut->data;
   for (uint32_t i = 2; i < nor; ++i) {
      for (uint32_t k = 1; k < i; ++k) {
         for (uint32_t l = 0; l < k; ++l) {
            p2[mape3(i, k, l)] = mape3(p1[i], p1[k], p1[l]);
         }
      }
   }
   permSave(permOut, fileNameOut);
}


/* ------------------------------------------------------------------
   ze4() - Antisymmetric fourth power
   ------------------------------------------------------------------ */

static void ze4()
{
   FEL f11, f12, f13, f14, f21, f22, f23, f24, f31, f32, f33, f34, f41, f42, f43, f44;
   FEL e, g12, g13, g14, g23, g24, g34, g123, g124, g134, g234;
   int i1, i2, i3, i4, j1, j2, j3, j4, jins;

   for (i1 = 0; i1 < nor - 3; ++i1) {
      MESSAGE(1, ("i1 = %d\n", i1));
      for (i2 = i1 + 1; i2 < nor - 2; ++i2) {
         MESSAGE(2, ("i2 = %d\n", i2));
         for (i3 = i2 + 1; i3 < nor - 1; ++i3) {
            MESSAGE(3, ("i3 = %d\n", i3));
            for (i4 = i3 + 1; i4 < nor; ++i4) {
               ffMulRow(rowOut, FF_ZERO, nocOut);
               jins = 0;
               for (j1 = 0; j1 < noc - 3; ++j1) {
                  f11 = ffExtract(rowIn[i1], j1);
                  f21 = ffExtract(rowIn[i2], j1);
                  f31 = ffExtract(rowIn[i3], j1);
                  f41 = ffExtract(rowIn[i4], j1);

                  for (j2 = j1 + 1; j2 < noc - 2; ++j2) {
                     f12 = ffExtract(rowIn[i1], j2);
                     f22 = ffExtract(rowIn[i2], j2);
                     f32 = ffExtract(rowIn[i3], j2);
                     f42 = ffExtract(rowIn[i4], j2);

                     g12 = ffSub(ffMul(f11, f22), ffMul(f21, f12));
                     g13 = ffSub(ffMul(f11, f32), ffMul(f31, f12));
                     g14 = ffSub(ffMul(f11, f42), ffMul(f41, f12));
                     g23 = ffSub(ffMul(f21, f32), ffMul(f31, f22));
                     g24 = ffSub(ffMul(f21, f42), ffMul(f41, f22));
                     g34 = ffSub(ffMul(f31, f42), ffMul(f41, f32));

                     for (j3 = j2 + 1; j3 < noc - 1; ++j3) {
                        f13 = ffExtract(rowIn[i1], j3);
                        f23 = ffExtract(rowIn[i2], j3);
                        f33 = ffExtract(rowIn[i3], j3);
                        f43 = ffExtract(rowIn[i4], j3);

                        g123 = ffMul(f13, g23);
                        g123 = ffSub(g123, ffMul(f23, g13));
                        g123 = ffAdd(g123, ffMul(f33, g12));
                        g124 = ffMul(f13, g24);
                        g124 = ffSub(g124, ffMul(f23, g14));
                        g124 = ffAdd(g124, ffMul(f43, g12));
                        g134 = ffMul(f13, g34);
                        g134 = ffSub(g134, ffMul(f33, g14));
                        g134 = ffAdd(g134, ffMul(f43, g13));
                        g234 = ffMul(f23, g34);
                        g234 = ffSub(g234, ffMul(f33, g24));
                        g234 = ffAdd(g234, ffMul(f43, g23));

                        for (j4 = j3 + 1; j4 < noc; ++j4) {
                           f14 = ffExtract(rowIn[i1], j4);
                           f24 = ffExtract(rowIn[i2], j4);
                           f34 = ffExtract(rowIn[i3], j4);
                           f44 = ffExtract(rowIn[i4], j4);

                           e = ffMul(f24, g134);
                           e = ffSub(e, ffMul(f14, g234));
                           e = ffAdd(e, ffMul(f44, g123));
                           e = ffSub(e, ffMul(f34, g124));
                           ffInsert(rowOut, jins, e);
                           ++jins;
                        }
                     }
                  }
               }
               mfWriteRows(fileOut, rowOut, 1, nocOut);
            }
         }
      }
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int init(int argc, char** argv)
{
   const char* arg3;

   /* Process command line options.
      ----------------------------- */
   App = appAlloc(&AppInfo, argc, argv);
   if (App == NULL) {
      return -1;
   }
   opt_G = appGetOption(App, "-G --gap");
   if (opt_G) {
      MtxMessageLevel = -100;
   }

   /* Process arguments.
      ------------------ */
   if (appGetArguments(App, 3, 3) < 0) {
      return -1;
   }
   fileNameInp = App->argV[1];
   fileNameOut = App->argV[2];
   arg3 = App->argV[0];
   if (!strcmp(arg3, "e2")) { mode = M_E2;} else if (!strcmp(arg3, "e3")) {
      mode = M_E3;
   } else if (!strcmp(arg3, "e4")) { mode = M_E4;} else if (!strcmp(arg3, "s2")) {
      mode = M_S2;
   } else {
      mtxAbort(MTX_HERE, "Unknown mode '%s'", arg3);
      return -1;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
   init(argc, argv);
   prepare();
   switch (mode) {
      case M_S2: objectType == MTX_TYPE_MATRIX ? zs2() : zs2p();
         break;
      case M_E2: objectType == MTX_TYPE_MATRIX ? ze2() : ze2p();
         break;
      case M_E3: objectType == MTX_TYPE_MATRIX ? ze3() : ze3p();
         break;
      case M_E4: ze4();
         break;
      default:
         mtxAbort(MTX_HERE, "Unknown mode %d", (int)mode);
         break;
   }
   if (fileOut != NULL) {
      mfClose(fileOut);
   }
   appFree(App);
   return 0;
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
