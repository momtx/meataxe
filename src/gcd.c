////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Seed vector generator
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>

/// @defgroup misc Miscellaneous
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculates the greatest common divisor of two integers.
/// If both numbers are positive, the result is positive. If one or both numbers are negative, the
/// result may be positive or negative. If one of the numbers is zero the function returns the
/// other number. If both numbers are zero, the result is zero.
/// @param a First number.
/// @param b Second number.
/// @return Greatest common divisor of @a a and @a b.

long gcd(long a,long b)
{
   if (a == 0) {
      return b;
   }
   if (b == 0) {
      return a;
   }
   while ((a %= b) != 0) {
      if ((b %= a) == 0) {
         return a;
      }
   }
   return b;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculates the least common multiple of two integers.
/// If both numbers are positive, the result is positive. If one or both numbers are negative,
/// the result may be positive or negative. If both numbers are zero, a division by zero occurs.
/// @param a First number.
/// @param b Second number.
/// @return Least common multiple of @a a and @a b.

long lcm(long a,long b)
{
   long ggt = gcd(a,b);
   if (ggt == 0) {
      return 0;
   }
   return (a / ggt) * b;
}


/// @}
