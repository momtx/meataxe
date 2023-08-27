////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Basic greased matrix functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


#define GMAT_MAGIC 0x52068001

/// @defgroup grmat Greased Matrices
/// @details
/// A greased matrix is a matrix over a finite field, which has been
/// optimized for fast row operations. The optimization ("grease") is
/// achieved by precomputing linear combinations of blocks of rows.
///
/// The number of rows per block, also called grease level, is
/// restricted to the range 1...16. Grease level 3, for example, means
/// that the rows of the matrix are divided in blocks of three rows, and
/// for each block, all linear combinations of the three rows are calculated
/// once. Multiplying a single vector by the matrix can then be carried out
/// with only n/3 row operations.
///
/// On the other hand, the greased matrix needs more memory. For grease level
/// 8 with GF(2), the memory needed is increased by a factor of 32.
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check a greased matrix.
/// This function checks if the argument |mat| is a pointer to a valid
/// greased matrix. If the matrix is o.k., the function returns 1.
/// Otherwise, an error is signalled and, if the error handler does not
/// terminate the program, the function returns 0.
/// @param mat Pointer to the matrix.
/// @return 1 if @a mat points to a valid greased matrix, 0 otherwise.

int GrMatIsValid(const GreasedMatrix_t *mat)
{
   if (mat == NULL) {
      mtxAbort(MTX_HERE,"NULL matrix");
      return 0;
   }
   if ((mat->Magic != GMAT_MAGIC) || (mat->Field < 2) || (mat->Nor < 0) ||
       (mat->Noc < 0)) {
      mtxAbort(MTX_HERE,"Invalid greased matrix (field=%d, nor=%d, noc=%d)",
                 mat->Field,mat->Nor,mat->Noc);
      return 0;
   }
   return 1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Free a greased matrix.
/// 0 on success, -1 on error.
/// This function frees a greased matrix, releasing all internally
/// allocated memory. Note that some data structures (the extraction
/// tables) are kept in a cache and are never freed until the process
/// terminates.
/// @param mat The matrix to be freed.
/// @return 0 on success, -1 on error.

int GrMatFree(GreasedMatrix_t *mat)
{
   if (!GrMatIsValid(mat)) {
      return -1;
   }
   if (mat->PrecalcData != NULL) {
      sysFree(mat->PrecalcData);
   }
   memset(mat,0,sizeof(Matrix_t));
   sysFree(mat);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Create a greased matrix.
/// This function creates a greased matrix from an existing normal matrix.
/// Basically, this means that linear combinations of the rows of @a mat
/// are calculated and stored in a buffer. The number of precalculated
/// rows depends on the field order and the grease level.
/// The original matrix is unchanged, and the caller is responsible for
/// deleting it. @a gr_nrows is the grease level, or block size, and must
/// be in the range 0...16. A grease level of 0 means that greasing is
/// switched off.
///
/// To destroy a greased matrix, use GrMatFree().
/// @param M The normal matrix.
/// @param gr_rows Grease level (number of rows per block).
/// @return Pointer to the greased matrix or 0 on error.

GreasedMatrix_t *GrMatAlloc(const Matrix_t *M, int gr_rows)
{
   long i,j,k,l;  /* counters */
   long rows;     /* number of rows written so far in current part */
   PTR p;         /* to write vectors */
   PTR q;         /* to read vectors from the original matrix */
   PTR r;         /* to read vectors from previous results */
   PTR bs;        /* start of current block of vectors */
   PTR v;         /* space for one vector */
   long nrvecs;   /* number of vectors in greased table */
   GreasedMatrix_t *res;   /* the result */

   res = ALLOC(GreasedMatrix_t);
   if (res == NULL) {
      return NULL;
   }
   res->Field = M->Field;
   res->Noc = M->Noc;
   res->Nor = M->Nor;
   res->GrRows = gr_rows;

   ffSetField(M->Field);

   /* special case of greasing switched off: */
   if (gr_rows == 0) {
      res->GrBlockSize = 1;
      res->PrecalcData = ffAlloc(M->Nor, M->Noc);
      res->NumVecs = M->Nor;
      memcpy(res->PrecalcData,M->Data,ffSize(M->Nor, M->Noc));
      res->ExtrTab = NULL;
      res->Magic = GMAT_MAGIC;
      return res;
   }

   nrvecs = 1;
   for (i = gr_rows; i > 0; i--) {
      nrvecs *= M->Field;   /* This calculates fl^gr_rows */
   }
   res->GrBlockSize = --nrvecs;   /* the null vector is not stored! */

   nrvecs = nrvecs * (M->Nor / gr_rows) + M->Nor % gr_rows;
   res->NumVecs = nrvecs;
   if ((res->PrecalcData = ffAlloc(nrvecs, M->Noc)) == NULL) {
      return NULL;
   }
   v = ffAlloc(1, M->Noc);

   /* Now we calculate all linear combinations necessary: */
   p = res->PrecalcData;
   q = M->Data;
   for (i = M->Nor / gr_rows; i > 0; i--) { /* all greasing block */
      bs = p;
      rows = 0;    /* this is 1 - 1 = fl^0 - 1 */

      /* now for all vectors in the block: */
      for (j = gr_rows; j > 0; j--) { /* all vectors in block */
         for (k = 1; k < M->Field; k++) { /* all field elements */
            ffCopyRow(v,q, M->Noc);
            ffMulRow(v,ffFromInt(k), M->Noc);
            ffCopyRow(p,v, M->Noc); /* copy the new multiple */
            ffStepPtr(&p, M->Noc);

            r = bs;  /* start from the beginning of the current block */
            for (l = rows; l > 0; l--) { /* for all vectors so far */
               ffCopyRow(p,r, M->Noc);
               ffStepPtr(&r, M->Noc);
               ffAddRow(p,v, M->Noc);
               ffStepPtr(&p, M->Noc);
            }
         }
         ffStepPtr(&q, M->Noc); /* take a new row of the original matrix */
         rows = (rows + 1) * M->Field - 1; /* the null vector is not there */
      }
   }

   for (i = M->Nor % gr_rows; i > 0; i--) { /* the rest of the vectors */
      ffCopyRow(p,q, M->Noc);
      ffStepPtr(&p, M->Noc);
      ffStepPtr(&q, M->Noc);
   }
   res->ExtrTab = GrGetExtractionTable(M->Field,gr_rows);
   sysFree(v);

   res->Magic = GMAT_MAGIC;
   return res;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
