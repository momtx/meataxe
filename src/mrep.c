////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix representations, core functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

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

static void validateGenerators(
      const struct MtxSourceLocation* where, int ngen, Matrix_t **gen)
{
   int i;

   if (ngen < 0) {
      mtxAbort(where,"ngen: %s",MTX_ERR_BADARG);
   }
   if ((ngen > 0) && (gen == NULL)) {
      mtxAbort(where,"gen == NULL: %s",MTX_ERR_BADARG);
   }
   for (i = 0; i < ngen; ++i) {
      matValidate(where, gen[i]);
      if (gen[i]->nor != gen[i]->noc) {
         mtxAbort(where,"gen[%i]: %s",i,MTX_ERR_NOTSQUARE);
      }
      if (i != 0) {
         if ((gen[i]->field != gen[0]->field) || (gen[i]->nor != gen[0]->nor)) {
            mtxAbort(where,"gen[0] and gen[%d]: %s",i,MTX_ERR_INCOMPAT);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Verifies that a matrix representation is valid and aborts the program if it is not valid.

void mrValidate(const struct MtxSourceLocation* where, const MatRep_t *rep)
{
   if (rep == NULL)
      mtxAbort(where,"NULL representation");
   if (rep->typeId != MTX_TYPE_MATREP) 
      mtxAbort(where,"Invalid matrix representation (type=0x%lx)",(unsigned long)rep->typeId);
   if (rep->NGen < 0)
      mtxAbort(where,"Invalid number of generators (%d)",rep->NGen);
   validateGenerators(MTX_HERE, rep->NGen,rep->Gen);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates a matrix representation.
///
/// The matrices in @a gen must all be square, over the same field, and of the same size.
///
/// @a flags may be zero or the special value MR_COPY_GENERATORS. In the latter case, a local copy
/// of the generators is made, and the matrices in @a gen can be safely destroyed. If @a flags is 0,
/// the representation becomes owner of the genertors, and the caller must not modify or free any
/// of them. The generators will be freed by @ref mrFree.

/// @param ngen Number of generators in @a gen.
/// @param gen List of generators. May be NULL if @a ngen is 0.
/// @param flags Optional flags (see description).
/// @return Pointer to the new matrix representation

// TODO: add a MR_FOREIGN_GENERATORS to prevent mrFree() from destroying the generators

MatRep_t *mrAlloc(int ngen, Matrix_t **gen, int flags)
{
   validateGenerators(MTX_HERE, ngen, gen);
   MatRep_t*rep = (MatRep_t*) mmAlloc(MTX_TYPE_MATREP, sizeof(MatRep_t));
   rep->Gen = NALLOC(Matrix_t *,ngen);
   rep->NGen = ngen;
   for (int i = 0; i < ngen; ++i) {
      rep->Gen[i] = (flags & MR_COPY_GENERATORS) ?  matDup(gen[i]) : gen[i];
   }
   return rep;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Delete a matrix representation.
///
/// Note: All generators (added by @ref mrAlloc or @ref mrAddGenerator) are destroyed, even if the
/// MR_COPY_GENERATORS was used.

void mrFree(MatRep_t *rep)
{
   mrValidate(MTX_HERE,rep);
   for (int i = 0; i < rep->NGen; ++i)
      matFree(rep->Gen[i]);
   memset(rep->Gen,0,sizeof(Matrix_t *) * rep->NGen);
   sysFree(rep->Gen);
   rep->Gen = NULL;
   mmFree(rep, MTX_TYPE_MATREP);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Add a generator to a matrix representation.
///
/// The matrix must be square. If there are already generators in the representation, the new
/// matrix must be over the same field and have the same number of rows.
///
/// @a flags may be zero or the special value MR_COPY_GENERATORS. See @ref mrAlloc.

void mrAddGenerator(MatRep_t *rep, Matrix_t *gen, int flags)
{
    mrValidate(MTX_HERE,rep);
    if (gen->nor != gen->noc)
	mtxAbort(MTX_HERE,"gen: %s",MTX_ERR_NOTSQUARE);
    if (rep->NGen > 0)
    {
	if (gen->field != rep->Gen[0]->field || gen->nor != rep->Gen[0]->nor)
	    mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
    }

    rep->Gen = NREALLOC(rep->Gen, Matrix_t *,rep->NGen + 1);
    rep->Gen[rep->NGen++] = (flags & MR_COPY_GENERATORS) ? matDup(gen) : gen;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Load a Matrix Representation.
/// This function creates a new matrix representation and reads the generators
/// from files. Each generator ist expected in a different file. The file name
/// is constructed by appending ".1", ".2" etc. to @a basename or, if @a basename
/// contains a "%d" placeholder, by replacing the "%d" with "1", "2", etc.
/// For example, the following lines
/// @code
/// m11 = mrLoad("m11",2);
/// m11 = mrLoad("m11.%d",2);
/// @endcode
/// are equivalent. In both cases, two matrices are read from "m11.2" and "m11.2",
/// repectively.
/// @param basename Base file name for generators.
/// @param ngen Number of generators.
/// @return Pointer to the representation.

MatRep_t *mrLoad(const char *basename, int ngen)
{
   // Make a copy of the basename and reserve extra bytes for the extension.
   char* fn = sysMalloc(strlen(basename) + 10);

   MatRep_t* mr = mrAlloc(0,NULL,0);

   // Read the generators
   int ext_format = strstr(basename,"%d") != NULL; // handle '%d' in <basename>
   for (int i = 0; i < ngen; ++i) {
      Matrix_t *gen;
      if (ext_format) {
         sprintf(fn,basename,i + 1);
      } else {
         sprintf(fn,"%s.%d",basename,i + 1);
      }
      gen = matLoad(fn);
      mrAddGenerator(mr,gen,0);
   }

   sysFree(fn);
   return mr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Transpose a Representation.
/// This function transposes a matrix representation. A new representation
/// is created, and the original representation is not changed.
/// @param rep Matrix representation.
/// @return Pointer to the transposed representation, or 0 on error.

MatRep_t *mrTransposed(const MatRep_t *rep)
{
   mrValidate(MTX_HERE,rep);

   Matrix_t **tr = NALLOC(Matrix_t *,rep->NGen);
   for (int i = 0; i < rep->NGen; ++i) {
      tr[i] = matTransposed(rep->Gen[i]);
   }

   MatRep_t *tr_rep = mrAlloc(rep->NGen,tr,0);

   sysFree(tr);
   return tr_rep;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Save a Matrix Representation.
/// This function saves all generators of a matrix representation.
/// Each generator ist written to different file. The file name
/// is constructed by appending ".1", ".2" etc. to @a basename or, if
/// @a basename contains a "%d" placeholder, by replacing the "%d"
/// with "1", "2", etc.
/// @param rep Pointer to the matrix representation.
/// @param basename Base file name for generators.
/// @return 0 on success, -1 on error.

int mrSave(const MatRep_t *rep, const char *basename)
{
   char *fn;
   int ext_format;          /* '%d' found in <basename> */
   int i;

   /* Make a copy of the basename an reserve extra bytes for the extension
      -------------------------------------------------------------------- */
   fn = sysMalloc(strlen(basename) + 10);
   if (fn == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate buffer");
      return -1;
   }

   /* Write the generators.
      --------------------- */
   ext_format = strstr(basename,"%d") != NULL;
   for (i = 0; i < rep->NGen; ++i) {
      if (ext_format) {
         sprintf(fn,basename,i + 1);
      } else {
         sprintf(fn,"%s.%d",basename,i + 1);
      }
      matSave(rep->Gen[i],fn);
   }

   /* Clean up.
      --------- */
   sysFree(fn);
   return i >= rep->NGen ? 0 : -1;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
