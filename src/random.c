////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Seed vector generator
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <limits.h>     /* For LONG_MAX */

////////////////////////////////////////////////////////////////////////////////////////////////////
//  Local data

MTX_DEFINE_FILE_INFO

static long int randtbl[32] = {
   3,
   -851904987, -43806228, -2029755270, 1390239686, -1912102820,
   -485608943,1969813258, -1590463333,-1944053249,455935928,508023712,
   -1714531963, 1800685987, -2015299881, 654595283, -1149023258,
   -1470005550, -1143256056, -1325577603, -1568001885, 1275120390,
   -607508183, -205999574, -1696891592, 1492211999, -1528267240,
   -952028296, -189082757, 362343714, 1424981831, 2039449641,
};

static long int *fptr = &randtbl[4];
static long int *rptr = &randtbl[1];
static long int *state = &randtbl[1];
static long int *end_ptr = &randtbl[sizeof(randtbl) / sizeof(randtbl[0])];

/// @addtogroup misc
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Initializes the random number generator.
/// @param x Seed value.

void MtxRandomInit(unsigned x)
{
   register int i;
   state[0] = x;
   for (i = 1; i < 31; ++i) {
      state[i] = (1103515145 * state[i - 1]) + 12345;
   }
   fptr = &state[3];
   rptr = &state[0];
   for (i = 0; i < 10 * 31; ++i) {
      MtxRandom();
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Get the next (pseudo-)random number
/// This function returns a nonnegative (pseudo) random number.
/// @return A nonnegative number.

long MtxRandom(void)
{
   long int i;

   *fptr += *rptr;
   i = (*fptr >> 1) & LONG_MAX;         /* Chuck least random bit.  */
   ++fptr;
   if (fptr >= end_ptr) {
      fptr = state;
      ++rptr;
   } else {
      ++rptr;
      if (rptr >= end_ptr) {
         rptr = state;
      }
   }
   return i;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Random number.
/// This function returns a (pseudo) random number in the range 0...@a max.
/// The argument must be greater than 0.
/// @return Pseudo random number.

int MtxRandomInt(int max)
{
   MTX_VERIFY(max > 0);
   return (int) (MtxRandom() % max);
}


/// @}