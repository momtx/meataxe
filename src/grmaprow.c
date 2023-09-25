////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Greased matrix multiplication
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>
#include <string.h>


#if MTX_ZZZ == 0	// greasing is only available for the standard kernel

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup grmat
/// @{


static long *extractednrs = NULL;
static long extrtablen = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void ExtractNrs(PTR v,GreasedMatrix_t *M)
{
/* Builds a table of numbers extracted from the vector v. Stores it into
   extractednrs (global variable, see above). This is malloced anew if
   not long enough! */
   unsigned char *p = v;
   register long *q;
   register int curtab = 0;
   register int i = M->nor / M->grRows;

   if (extrtablen < i + 16) { /* some reserve! */
      if (extractednrs) { sysFree(extractednrs); }
      extractednrs = NALLOC(long,i + 16);
      if (extractednrs == NULL) {
         mtxAbort(MTX_HERE,"Not enough memory for extraction table");
         return;
      }
      extrtablen = i + 16;
   }
   q = extractednrs;

   {
      register int y = 0;
      register long *x;

      while (i > 0) {
         /* tabnr = 0; */
         x = M->extrTab->tabs[curtab][*p++];
         {
            register int k = M->extrTab->nrvals[curtab];
            while (k--) {
               *q++ = y + *x++;
               y = 0;
               i--;
            }
         }
         if (++curtab < M->extrTab->nrtabs) {
            y += *x; /* remember last part for "carry" */
         } else {
            curtab = 0;
         }
      }
   }
}


static void CleverMapRow(PTR v,GreasedMatrix_t *M,PTR w)
{
/* Calculates the matrix product of the vector `v' with the matrix `M'. The
   result is stored in `w'. The length of `v' must coincide with the number
   of rows of `M'. The result `w' is as long as the length of a row of
   `M'. This is called from `GreaseMapRow' in the case where we have
   an extraction table. */
   register long *nr;
   register long i;
   long curcol;
   PTR p = M->precalcData;

   ExtractNrs(v,M);    /* extracts numbers into global table */
   ffMulRow(w,FF_ZERO, M->noc);   /* Clear result */
   curcol = M->nor / M->grRows;   /* only used later and for i init */
   nr = extractednrs;
   for (i = curcol; i > 0; i--) {
      /* for all greasing blocks */
      if (*nr) { /* not the zero vector */
         ffAddRow(w,ffGetPtr(p,*nr - 1, M->noc), M->noc);
      }
      p = ffGetPtr(p,M->grBlockSize, M->noc);
      nr++;
   }
   curcol = curcol * M->grRows;
   for (i = M->nor % M->grRows; i > 0; i--) {
      /* do the rest */
      ffAddMulRow(w,p,ffExtract(v,curcol++),M->noc);
      ffStepPtr(&p, M->noc);
   }
}


/// Multiply a vector by a greased matrix.
/// This function calculates the matrix product of the vector @a v with
/// the matrix @a M and write the result to @a w. The length of @a v must
/// coincide with the number of rows of @a M. The result @a w is as long
/// as the length of a row of @a M.
/// Unlike ffMapRow(), this function sets field and row length correctly!
/// @param v The vector.
/// @param M Pointer to the matrix.
/// @param w The result, vM.

void GrMapRow(PTR v,GreasedMatrix_t *M, PTR w)
{
   int i,j;      /* Counters */
   int curcol;   /* Current column of `v' */
   int nr;       /* number of vector in the lookup table */
   PTR p;         /* Pointer into lookup table */

   if (!GrMatIsValid(M)) {
      mtxAbort(MTX_HERE,"GrMapRow(): Invalid argument(s)");
      return;
   }
   ffSetField(M->field);

   /* Handle some special cases.
      -------------------------- */
   if (M->extrTab) {            /* We use our clever extraction */
      CleverMapRow(v,M,w);
      return;
   }
   if (M->grRows == 0) {        /* Greasig is switched off */
      ffMapRow(v,M->precalcData,M->nor,M->noc,w);
      return;
   }

   ffMulRow(w,FF_ZERO, M->noc);   /* Clear result */
   curcol = 0;                  /* start with the first column of `v' */
   p = M->precalcData;          /* start at the beginning of the table */

   if (M->field == 2) {
      register signed char *q = (signed char *) v;
      register signed char buf = *q++;
      register long wmask;     /* used to alter nr */
      for (i = M->nor / M->grRows; i > 0; i--) {  /* for all greasing blocks */
         /* Calculate the number for the first lookup: */
         nr = 0;
         wmask = 1;
         for (j = M->grRows - 1; j >= 0; j--) {
            if (buf < 0) {    /* a bit is set */
               nr |= wmask;
            }
            wmask <<= 1;
            if ((++curcol & 0x7) == 0) {
               buf = *q++;
            } else {
               buf <<= 1;
            }
         }
         if (nr) {     /* not the zero vector */
            p = ffGetPtr(p,nr - 1,M->noc);   /* now p points to the right linear combination */
            ffAddRow(w,p, M->noc);
            p = ffGetPtr(p,M->grBlockSize - nr + 1, M->noc);
         } else {
            p = ffGetPtr(p,M->grBlockSize, M->noc);
         }
      }
      for (i = M->nor % M->grRows; i > 0; i--) {
         /* do the rest */
         if (buf < 0) { ffAddRow(w,p, M->noc); }
         ffStepPtr(&p, M->noc);
         if ((++curcol & 0x7) == 0) {
            buf = *q++;
         } else {
            buf <<= 1;
         }
      }
   } else {
      /* not GF(2) */
      for (i = M->nor / M->grRows; i > 0; i--) {  /* for all greasing blocks */
         /* Calculate the number for the first lookup: */
         nr = 0;
         for (j = M->grRows - 1; j >= 0; j--) { /* add j to curcol */
            nr = nr * M->field + ffToInt(ffExtract(v,curcol + j));
         }
         if (nr) {     /* not the zero vector */
            p = ffGetPtr(p,nr - 1, M->noc);   /* now p points to the right linear combination */
            ffAddRow(w,p, M->noc);
            p = ffGetPtr(p,M->grBlockSize - nr + 1, M->noc);
         } else {
            p = ffGetPtr(p,M->grBlockSize, M->noc);
         }
         curcol += M->grRows;
      }
      for (i = M->nor % M->grRows; i > 0; i--) {  /* do the rest */
         ffAddMulRow(w,p,ffExtract(v,curcol++),M->noc);
         ffStepPtr(&p, M->noc);
      }
   }
}

#else
// greasing is only available for the standard kernel

void GrMapRow(PTR v,GreasedMatrix_t *M, PTR w)
{
    mtxAbort(MTX_HERE,"Greasing is not yet supported for ZZZ=1");
    abort();
}

#endif

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
