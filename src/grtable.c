////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculation of extraction tables for greasing
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup grmat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Builds the tables used for the extraction of grease bits.
//  It can build tables for the following field/grease combinations: field: 2-256, grrows: 1-16

static GrExtractionTable_t *BuildExtractionTable(int fl,int grrows)
{
   int MPB;                     /* marks per byte */
   GrExtractionTable_t *t;      /* the result */
   unsigned char c[4];          /* buffer for calculations */
   long flpowtab[17];           /* powers of fl */
   long *w;                     /* to write into tables */
   int i,j,k,l;                 /* counters */
   int nrvals;
   int restbits;
   static const char no_mem[] = "Not enough memory for extraction table";

   if ((t = ALLOC(GrExtractionTable_t)) == NULL) {
      mtxAbort(MTX_HERE,no_mem);
      return NULL;
   }

   /* Note that sizeof(long)*2*8*5*3 will be divisable by MPB*sizeof(long)*2,
      so a vector of this length will completely use all zsize(..) bytes.
      We need the *2 because of possible MMX modifications. */
   ffSetField(fl);
   int noc = sizeof(long) * 2 * 8 * 5 * 3;
   ffSetNoc(noc);
   MPB = noc / ffRowSize(noc);

   /* Allocate memory.
      ---------------- */
   l = lcm(MPB,grrows);         /* Number of byte types (=number of tables) */
   t->nrtabs = l;
   t->tabs = NALLOC(long **,l);
   t->nrvals = NALLOC(int,l);
   if ((t->tabs == NULL) || (t->nrvals == NULL)) {
      mtxAbort(MTX_HERE,no_mem);
      return NULL;
   }
   /* calculate flpowtab: */
   for (l = 0, j = 1; l <= 16; l++, j *= fl) {
      flpowtab[l] = j;
   }

   for (l = 0; l < t->nrtabs; l++) {
      if ((t->tabs[l] = NALLOC(long*,256)) == NULL) {
         mtxAbort(MTX_HERE,no_mem);
         return NULL;
      }

      /* Calculate the number of values completed in this byte: */
      nrvals = 0;     /* we count the number of values */
      if ((l * MPB) % grrows > 0) {
         restbits = grrows - (l * MPB) % grrows;
         if (restbits <= MPB) {    /* one value to complete: */
            nrvals++;
         }
      } else {
         restbits = 0;
      }
      nrvals += (MPB - restbits) / grrows;
      t->nrvals[l] = nrvals;

      /* Now we go through all possibilities of this byte: */
      for (i = 0; i < flpowtab[MPB]; i++) {
         *c = 0;
         for (j = 1,k = i; j <= MPB; j++) {
            ffInsert((PTR)c,j - 1,ffFromInt(k % fl));
            k /= fl;
         }

         if ((l < 0) || (l > t->nrtabs)) {
            mtxAbort(MTX_HERE,"Invalid table number %d",l);
         }
         if ((t->tabs[l][*c] = NALLOC(long,nrvals + 1)) == NULL) {
            mtxAbort(MTX_HERE,no_mem);
            return NULL;
         }

         /* now distribute the different values in i to the table: */
         w = t->tabs[l][*c];
         k = i;
         /* first (perhaps) a value that already begun in the previous byte: */
         if (restbits > 0) { /* one value to complete: */
            if (w - t->tabs[l][*c] > nrvals + 1) {
               mtxAbort(MTX_HERE,"Table overflow");
            }
            if (restbits <= MPB) {
               /* is completed within this byte! */
               *w++ = flpowtab[(l * MPB) % grrows] * (k % flpowtab[restbits]);
               k /= flpowtab[restbits];
            } else {
               /* has already started and does not end within this byte: */
               *w++ = flpowtab[(l * MPB) % grrows] * k;
               continue;   /* to avoid the code for the last value */
            }
         }
         /* Now all the values fully within the byte: */
         for (j = (MPB - restbits) / grrows; j > 0; j--) {
            if (w - t->tabs[l][*c] > nrvals + 1) {
               mtxAbort(MTX_HERE,"Table overflow");
            }
            *w++ = k % flpowtab[grrows];
            k /= flpowtab[grrows];
         }
         /* Now (perhaps) a value that is only begun: */
         if (((l + 1) * MPB) % grrows > 0) { /* one additional part-value */
            if (w - t->tabs[l][*c] > nrvals + 1) {
               mtxAbort(MTX_HERE,"Table overflow");
            }
            *w++ = k;
         } else {
            if (w - t->tabs[l][*c] > nrvals + 1) {
               mtxAbort(MTX_HERE,"Table overflow");
            }
            *w++ = 0;
         }
      }
   }
   return t;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculate extraction table for grease.
/// This function calculates the extraction table for greased matrix
/// operations for a particular combination of field order and grease level.
/// The grease level must be in the range 1...16.
/// To avoid frequent table recalculations, tables are stored in a cache.
/// @param fl Field order.
/// @param grrows Grease level, number of rows per block.
/// @return Pointer to the extraction table, or NULL on error.

const GrExtractionTable_t *GrGetExtractionTable(int fl,int grrows)
{
   static GrExtractionTable_t *cache[257][17] = {{0}};

   if ((fl < 2) || (fl > 256)) {
      mtxAbort(MTX_HERE,"Invalid field order %d",fl);
      return NULL;
   }
   if ((grrows < 1) || (grrows > 16)) {
      mtxAbort(MTX_HERE,"Invalid grease level %d",grrows);
      return NULL;
   }

   if (cache[fl][grrows] == NULL) {
      cache[fl][grrows] = BuildExtractionTable(fl,grrows);
   }

   return cache[fl][grrows];
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
