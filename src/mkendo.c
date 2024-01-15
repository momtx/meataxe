////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Find endomorphisms
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <inttypes.h>


/// @defgroup endo Endomorphisms
/// @{

static Matrix_t *MakeEndo(const MatRep_t *rep, const Matrix_t *sb1, 
    const Matrix_t *vec)
{
    Matrix_t *sb2;
    Matrix_t *endo;

    // Make standard basis from <vec>
    sb2 = spinupStandardBasis(NULL, vec,rep,SF_FIRST);
    MTX_ASSERT(sb2 != NULL && sb2->nor == sb2->noc);
    //Matrix_t* sb2OLD = SpinUp(vec,rep,SF_FIRST|SF_CYCLIC|SF_STD,NULL,NULL);
    //MTX_ASSERT(matCompare(sb2, sb2OLD) == 0);
    //fprintf(stderr, "spinup %"PRIu32"x%"PRIu32"\n", sb2->nor, sb2->noc);
    //matFree(sb2OLD);

    // The linear mapping that maps <sb1> on <Matrix_t> is the endomorphism we are looking for!
    endo = matInverse(sb1);
    matMul(endo,sb2);
    matFree(sb2);
    return endo;
}



////////////////////////////////////////////////////////////////////////////////////////////////////

/// Calculates the endomorphism ring of an irreducible module.
///
/// The endomorphism ring E of an irreducible module has dimension equal to the degree of the
/// splitting field extension, d = [F':F].
///
/// On success, basis of E (a set of matrices) is stored into @a endo.
/// If the function fails, no matrices are stored in @a endo. 
/// @param rep
///   Pointer to an irreducible matrix representation.
/// @param nsp
///   Idword kernel. The number of rows must be equal to d.
///   extension.
/// @param endo
///   Result buffer. After successful execution, @p endo contains a basis of the 
///   endomorphism ring. The buffer size must be equal to d.
/// @return
///   0 on success, -1 on error.
///   If the return value is 0, @p endo contains d matrices which must be destroyed by the caller.
///   If an error occurred, the result buffer may have been modified, but it contains no valid
///   pointers.

int makeEndomorphisms(const MatRep_t* rep, const Matrix_t* nsp, Matrix_t* endo[])
{
   MTX_ASSERT(nsp->nor > 0);
   MTX_ASSERT(rep->NGen > 0);

   // Take the first vector from <nsp> and make the standard basis.
   Matrix_t* sb1 = spinupStandardBasis(NULL, nsp, rep, SF_FIRST);
   MTX_ASSERT(sb1 != NULL && sb1->nor == sb1->noc);
   //fprintf(stderr, "MKENDO %" PRIu32 "x%" PRIu32 "\n", nsp->nor, nsp->noc);
   //Matrix_t* sb1OLD = SpinUp(nsp, rep, SF_FIRST | SF_CYCLIC | SF_STD, NULL, NULL);
   //MTX_ASSERT(matCompare(sb1, sb1OLD) == 0);
   //matFree(sb1OLD);

   // Take the identity as the first basis element for E
   endo[0] = matId(rep->Gen[0]->field, rep->Gen[0]->nor);
   uint32_t nendo = 1;

   // For each of the remaining vectors v_2,..v_d in <nsp>, construct the
   // endomorphism that maps v_1 to v_j.
   while (nendo < nsp->nor) {
      Matrix_t* vec = matDupRows(nsp, nendo, 1);
      endo[nendo] = MakeEndo(rep, sb1, vec);
      matFree(vec);
      if (endo[nendo] == NULL) {        /* Error */
         break;
      }
      ++nendo;
   }
   matFree(sb1);

   // Clean up after error.
   if (nendo < nsp->nor) {
      while (nendo > 0) {
         matFree(endo[--nendo]);
      }
   }

   return nendo == nsp->nor ? 0 : -1;
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
