////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Greased matrix multiplication
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

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
   register int i = M->Nor / M->GrRows;

   if (extrtablen < i + 16) { /* some reserve! */
      if (extractednrs) { SysFree(extractednrs); }
      extractednrs = NALLOC(long,i + 16);
      if (extractednrs == NULL) {
         MTX_ERROR("Not enough memory for extraction table");
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
         x = M->ExtrTab->tabs[curtab][*p++];
         {
            register int k = M->ExtrTab->nrvals[curtab];
            while (k--) {
               *q++ = y + *x++;
               y = 0;
               i--;
            }
         }
         if (++curtab < M->ExtrTab->nrtabs) {
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
   PTR p = M->PrecalcData;

   ExtractNrs(v,M);    /* extracts numbers into global table */
   FfMulRow(w,FF_ZERO);   /* Clear result */
   curcol = M->Nor / M->GrRows;   /* only used later and for i init */
   nr = extractednrs;
   for (i = curcol; i > 0; i--) {
      /* for all greasing blocks */
      if (*nr) { /* not the zero vector */
         FfAddRow(w,FfGetPtr(p,*nr - 1));
      }
      p = FfGetPtr(p,M->GrBlockSize);
      nr++;
   }
   curcol = curcol * M->GrRows;
   for (i = M->Nor % M->GrRows; i > 0; i--) {
      /* do the rest */
      FfAddMulRow(w,p,FfExtract(v,curcol++));
      FfStepPtr(&p);
   }
}


/// Multiply a vector by a greased matrix.
/// This function calculates the matrix product of the vector @a v with
/// the matrix @a M and write the result to @a w. The length of @a v must
/// coincide with the number of rows of @a M. The result @a w is as long
/// as the length of a row of @a M.
/// Unlike FfMapRow(), this function sets field and row length correctly!
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
      MTX_ERROR("GrMapRow(): Invalid argument(s)");
      return;
   }
   FfSetField(M->Field);
   FfSetNoc(M->Noc);

   /* Handle some special cases.
      -------------------------- */
   if (M->ExtrTab) {            /* We use our clever extraction */
      CleverMapRow(v,M,w);
      return;
   }
   if (M->GrRows == 0) {        /* Greasig is switched off */
      FfMapRow(v,M->PrecalcData,M->Nor,w);
      return;
   }

   FfMulRow(w,FF_ZERO);   /* Clear result */
   curcol = 0;                  /* start with the first column of `v' */
   p = M->PrecalcData;          /* start at the beginning of the table */

   if (M->Field == 2) {
      register signed char *q = (signed char *) v;
      register signed char buf = *q++;
      register long wmask;     /* used to alter nr */
      for (i = M->Nor / M->GrRows; i > 0; i--) {  /* for all greasing blocks */
         /* Calculate the number for the first lookup: */
         nr = 0;
         wmask = 1;
         for (j = M->GrRows - 1; j >= 0; j--) {
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
            p = FfGetPtr(p,nr - 1);   /* now p points to the right linear combination */
#if 1
            FfAddRow(w,p);
#else
            {
               register long *a = (long *) w;
               register long *b = (long *) p;
               register long k = zsize(1) / sizeof(long);
               for (; k; k--) {
                  *a++ ^= *b++;
               }
            }
#endif
            p = FfGetPtr(p,M->GrBlockSize - nr + 1);
         } else {
            p = FfGetPtr(p,M->GrBlockSize);
         }
      }
      for (i = M->Nor % M->GrRows; i > 0; i--) {
         /* do the rest */
         if (buf < 0) { FfAddRow(w,p); }
         FfStepPtr(&p);
         if ((++curcol & 0x7) == 0) {
            buf = *q++;
         } else {
            buf <<= 1;
         }
      }
   } else {
      /* not GF(2) */
      for (i = M->Nor / M->GrRows; i > 0; i--) {  /* for all greasing blocks */
         /* Calculate the number for the first lookup: */
         nr = 0;
         for (j = M->GrRows - 1; j >= 0; j--) { /* add j to curcol */
            nr = nr * M->Field + FfToInt(FfExtract(v,curcol + j));
         }
         if (nr) {     /* not the zero vector */
            p = FfGetPtr(p,nr - 1);   /* now p points to the right linear combination */
            FfAddRow(w,p);
            p = FfGetPtr(p,M->GrBlockSize - nr + 1);
         } else {
            p = FfGetPtr(p,M->GrBlockSize);
         }
         curcol += M->GrRows;
      }
      for (i = M->Nor % M->GrRows; i > 0; i--) {  /* do the rest */
         FfAddMulRow(w,p,FfExtract(v,curcol++));
         FfStepPtr(&p);
      }
   }
}


/// @}