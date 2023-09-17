////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix representations, core functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
///   Local data


#define MR_MAGIC 0x1bb50442

/// @defgroup mrep Matrix Representations
/// @{
/// @class MatRep_t
/// A matrix representation.
/// This data structure is basically a collection of N matrices, representing the generators
/// of the algebra. A matrix representation is created with mrAlloc(). If the generators are
/// available, they can be passed to mrAlloc(). Generators can also be added to a matrix
/// representation with mrAddGenerator(). In both cases the user can choose if a copy of the
/// matrices or only a reference to the matrices is stored in the representation. In either
/// case, deleting the representation with mrFree() will also delete the generators.
/// For example, after the following code has been executed:
/// @code
/// Matrix_t *mat = matLoad("matrix");
/// Rep = mrAlloc(1,&mat,0);
/// mrFree(Rep);                               // *** invalidates mat
/// @endcode
/// the @c mat pointer is no longer valid, because the matrix was deleted by mrFree().
/// However, after
/// @code
/// Matrix_t *mat = matLoad("matrix");
/// Rep = mrAlloc(1,&mat,MR_COPY_GENERATORS);
/// mrFree(Rep);                               // *** mat remains valid
/// @endcode
/// @c mat is still valid because a copy of the matrix was created.

////////////////////////////////////////////////////////////////////////////////////////////////////

static int GensAreValid(int ngen, Matrix_t **gen)
{
   int i;

   if (ngen < 0) {
      mtxAbort(MTX_HERE,"ngen: %s",MTX_ERR_BADARG);
      return 0;
   }
   if ((ngen > 0) && (gen == NULL)) {
      mtxAbort(MTX_HERE,"gen == NULL: %s",MTX_ERR_BADARG);
      return 0;
   }
   for (i = 0; i < ngen; ++i) {
      matValidate(MTX_HERE, gen[i]);
      if (gen[i]->Nor != gen[i]->Noc) {
         mtxAbort(MTX_HERE,"gen[%i]: %s",i,MTX_ERR_NOTSQUARE);
         return 0;
      }
      if (i != 0) {
         if ((gen[i]->Field != gen[0]->Field) || (gen[i]->Nor != gen[0]->Nor)) {
            mtxAbort(MTX_HERE,"gen[0] and gen[%d]: %s",i,MTX_ERR_INCOMPAT);
            return 0;
         }
      }
   }
   return 1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check a matrix representation.
/// This function checks if the argument is a pointer to a valid
/// matrix representation. If the representation is o.k., the function
/// returns 1. Otherwise, an error is signaled and, if the error handler
/// does not terminate the program, the function returns 0.
/// @param rep Pointer to the matrix representation.
/// @return 1 if @a rep points to a valid matrix representation, 0 otherwise.

int mrIsValid(const MatRep_t *rep)
{
   if (rep == NULL) {
      mtxAbort(MTX_HERE,"NULL representation");
      return 0;
   }
   if (rep->Magic != MR_MAGIC) {
      mtxAbort(MTX_HERE,"Invalid matrix representation (magic=%d)",(int)rep->Magic);
      return 0;
   }
   if (!GensAreValid(rep->NGen,rep->Gen)) {
      mtxAbort(MTX_HERE,"Invalid generators");
      return 0;
   }
   return 1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Create a matrix representation.
/// The matrices in @a gen must all be square, over the same field, and of the same size.
/// @a flags may be zero or the special value MR_COPY_GENERATORS. In the latter case,
/// a local copy of the generators is made, and the matrices in @a gen can
/// be safely destroyed. If @a flags is 0, only references to the matrices
/// are stored in the MatRep_t structure. Consequently, the application
/// must not modify or destroy the matrices after calling mrAlloc(). They
/// will be destroyed automatically when mrFree() is called to destroy the
/// representation.
/// @param ngen Number of generators in @a gen.
/// @param gen List of generators. May be NULL if @a ngen is 0.
/// @param flags Optional flags (see description).
/// @return Pointer to the new matrix representation

// TODO: add a MR_FOREIGN_GENERATORS to prevent mrFree() from destroying the generators

MatRep_t *mrAlloc(int ngen, Matrix_t **gen, int flags)
{
   MatRep_t *rep;
   int i;

   if (!GensAreValid(ngen,gen)) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_BADARG);
      return NULL;
   }

   // Allocate a new MatRep_t structure
   rep = ALLOC(MatRep_t);
   if (rep == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate MatRep_t structure");
      return NULL;
   }
   memset(rep,0,sizeof(MatRep_t));
   rep->Gen = NALLOC(Matrix_t *,ngen);
   if (rep->Gen == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate generator list");
      sysFree(rep);
      return NULL;
   }

   // Copy generators
   rep->NGen = ngen;
   for (i = 0; i < ngen; ++i) {
      if (flags & MR_COPY_GENERATORS) {
         rep->Gen[i] = matDup(gen[i]);
         if (rep->Gen[i] == NULL) {
            mtxAbort(MTX_HERE,"Cannot copy generator");
            while (--i >= 0) {
               matFree(rep->Gen[i]);
            }
            sysFree(rep->Gen);
            sysFree(rep);
            return NULL;
         }
      } else {
         rep->Gen[i] = gen[i];
      }
   }

   rep->Magic = MR_MAGIC;
   return rep;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Delete a matrix representation.
/// This function frees a matrix representation which has beed created by
/// mrAlloc(). This implies freeing the internal data buffers as well as
/// the MatRep_t structure itself. Note: Even if the representation
/// was created without MR_COPY_GENERATORS, the matrices that were passed
/// to mrAlloc() are now destroyed. The same applies to matrices added to
/// the representation with mrAddGenerator().
/// @param rep Pointer to the matrix representation.
/// @return 0 on success, -1 on error.

int mrFree(MatRep_t *rep)
{
   int i;
   if (!mrIsValid(rep)) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_BADARG);
      return -1;
   }
   for (i = 0; i < rep->NGen; ++i) {
      matFree(rep->Gen[i]);
   }
   memset(rep->Gen,0,sizeof(Matrix_t *) * rep->NGen);
   sysFree(rep->Gen);
   memset(rep,0,sizeof(MatRep_t));
   sysFree(rep);
   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
