////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Hash function
// Based on lookup3.c by Bob Jenkins. Public Domain
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#define hashsize(n) ((uint32_t)1 << (n))
#define hashmask(n) (hashsize(n) - 1)
#define rot(x, k) (((x) << (k)) | ((x) >> (32 - (k))))

#define mix(a, b, c) \
   { \
      a -= c;  a ^= rot(c, 4);  c += b; \
      b -= a;  b ^= rot(a, 6);  a += c; \
      c -= b;  c ^= rot(b, 8);  b += a; \
      a -= c;  a ^= rot(c, 16);  c += b; \
      b -= a;  b ^= rot(a, 19);  a += c; \
      c -= b;  c ^= rot(b, 4);  b += a; \
   }

#define final(a, b, c) \
   { \
      c ^= b; c -= rot(b, 14); \
      a ^= c; a -= rot(c, 11); \
      b ^= a; b -= rot(a, 25); \
      c ^= b; c -= rot(b, 16); \
      a ^= c; a -= rot(c, 4);  \
      b ^= a; b -= rot(a, 14); \
      c ^= b; c -= rot(b, 24); \
   }


/// Returns two 32-bit hash values.
/// Both «pb» and «pc» are in-out parameters. The value passed in is used as seed. Each output
/// depends on both seeds.
/// *pc is better mixed than *pb, so use *pc first.
/// If you want a 64-bit value do something like "*pc + (((uint64_t)*pb)<<32).

void hashLittle2(const void *key, size_t length, uint32_t *pc, uint32_t *pb)
{
   // Set up the internal state
   uint32_t a, b, c; // internal state
   a = b = c = 0xdeadbeef + ((uint32_t)length) + *pc;
   c += *pb;

   {
      const uint8_t *k = (const uint8_t *)key;

      /*--------------- all but the last block: affect some 32 bits of (a,b,c) */
      while (length > 12)
      {
         a += k[0];
         a += ((uint32_t)k[1]) << 8;
         a += ((uint32_t)k[2]) << 16;
         a += ((uint32_t)k[3]) << 24;
         b += k[4];
         b += ((uint32_t)k[5]) << 8;
         b += ((uint32_t)k[6]) << 16;
         b += ((uint32_t)k[7]) << 24;
         c += k[8];
         c += ((uint32_t)k[9]) << 8;
         c += ((uint32_t)k[10]) << 16;
         c += ((uint32_t)k[11]) << 24;
         mix(a, b, c);
         length -= 12;
         k += 12;
      }

      // last block: affect all 32 bits of (c)
      switch (length)
      {
      case 12: c += ((uint32_t)k[11]) << 24; // fall-through
      case 11: c += ((uint32_t)k[10]) << 16; // fall-through
      case 10: c += ((uint32_t)k[9]) << 8; // fall-through
      case 9: c += k[8]; // fall-through
      case 8: b += ((uint32_t)k[7]) << 24; // fall-through
      case 7: b += ((uint32_t)k[6]) << 16; // fall-through
      case 6: b += ((uint32_t)k[5]) << 8; // fall-through
      case 5: b += k[4]; // fall-through
      case 4: a += ((uint32_t)k[3]) << 24; // fall-through
      case 3: a += ((uint32_t)k[2]) << 16; // fall-through
      case 2: a += ((uint32_t)k[1]) << 8; // fall-through
      case 1: a += k[0];
         break;
      case 0:
         *pc = c; *pb = b;
         return; // zero length strings require no mixing
      }
   }

   final(a, b, c);
   *pc = c;
   *pb = b;
}

// vim:et:sw=3:cin:fileencoding=utf8:fileformat=unix
