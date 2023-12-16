////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Find peak words and condense
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <string.h>

/// @addtogroup gap
/// @{

static uint32_t felToGapQ = 0;
static char** felToGapTable = NULL;
static char* felToGapBuffer = NULL;
static char primeFieldGenerator[15];

#define FEL_TO_INDEX(a) ((uint32_t) (a) >= felToGapQ ? felToGapQ-1 : (a))

static void rebuildTable()
{
    const size_t maxSize = ffOrder * 15;  // based on "Z(65536)^65534" (longest for ZZZ=1)
    felToGapBuffer = sysRealloc(felToGapBuffer, maxSize);
    char *wp = felToGapBuffer;
    char * const end = felToGapBuffer + maxSize;
    felToGapTable = NALLOC(char*, ffOrder);
    memset(felToGapTable, 0, ffOrder * sizeof(char*));
    felToGapQ = ffOrder;

    snprintf(primeFieldGenerator,sizeof(primeFieldGenerator),"Z(%u)", (unsigned) ffChar);
    if (ffOrder == ffChar) {
	// Prime field: "0", "1", "2", ... "q-1"
	FEL a = FF_ZERO;
	for (uint32_t k = 0; k < ffOrder; ++k) {
	    felToGapTable[FEL_TO_INDEX(a)] = wp;
	    wp += snprintf(wp, end - wp, "Z(%u)*%u", (unsigned)ffChar, (unsigned) k) + 1;
            MTX_ASSERT(wp <= end);
	    a = ffAdd(a, ffGen);
	}
	MTX_ASSERT(a == FF_ZERO);
    } else {
	// F(p^n), n>1
	// zero element: "0*Z(p)"
        felToGapTable[FEL_TO_INDEX(FF_ZERO)] = wp;
        wp += snprintf(wp, end - wp, "Z(%u)*0", (unsigned) ffOrder) + 1;
	// invertible elements: "Z(p)^k"
	FEL a = FF_ONE;
	for (uint32_t k = 0; k < ffOrder - 1; ++k) {
	    const uint32_t idx = FEL_TO_INDEX(a);
	    MTX_ASSERT(idx < ffOrder);
            MTX_ASSERT(felToGapTable[idx] == NULL);
            felToGapTable[idx] = wp;
            wp += snprintf(wp, end - wp, "Z(%u)^%u", (unsigned) ffOrder, (unsigned) k) + 1;
	    a = ffMul(a, ffGen);
	}
	MTX_ASSERT(a == FF_ONE);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns a GAP representation of a field element.
///
/// For prime fields, the returned string has the form "Z(p)*k" with 0≤k<p, where "Z(p)*1"
/// corresponds to the generator, @ref ffGen. Note that "Z(p)*0" is the zero element, but
/// "Z(p)*1" is not the unit element, except for p=2.
///
/// For fields of order q=p<sup>n</sup>, the zero element is represented as "Z(q)*0",
/// and nonzero elements are represented as "Z(q)^k", where "0≤k≤q-1". For example, if q = 25,
/// gapFelToString(FF_ZERO) returns "0*Z(25)", and gapFelToString(FF_ONE) returns "Z(25)^0".

const char* gapFelToString(FEL a)
{
   if (ffOrder != felToGapQ)
       rebuildTable();
   #if MTX_ZZZ == 0 || MTX_ZZZ == 1
   // ZZZ = 0: Numeric range of FEL is {0, 1, ... , q-1}
   // ZZZ = 1: Numeric range of FEL is {0, 1, ... , q-2, 0xFFFF}
   return ((uint32_t) a >= felToGapQ) ? felToGapTable[felToGapQ - 1] : felToGapTable[a];
   #else
      #error "Unsupported MTX_ZZZ"
   #endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the variable part of gapFelToString(). For example, if gapFelToString() would return
/// "Z(19)*7", gapFelToString1() returns "7".
/// The function fails and aborts the program if the current field order is not prime.
/// See also @ref gapFelToString2.

const char* gapFelToString1(FEL a)
{
   const char *s = strchr(gapFelToString(a), '*');
   if (s == NULL)
      mtxAbort(MTX_HERE, "%s(): argument 0x%x is not in GF(%u)",
           __func__, (unsigned)a, (unsigned) ffChar);
   return s + 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the fixed part of gapFelToString(). For example, if gapFelToString() would return
/// "Z(19)*7", gapFelToString2() returns "Z(19)".
/// The function fails and aborts the program if the current field order is not prime.
/// See also @ref gapFelToString2.

const char* gapFelToString2()
{
   if (ffOrder != felToGapQ)
      rebuildTable();
   if (ffOrder != ffChar)
      mtxAbort(MTX_HERE, "%s(): current field order %u is not prime", __func__, (unsigned)ffOrder);
   return primeFieldGenerator;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void gapFormatFel(StrBuffer_t* sb, FEL a)
{
    sbAppend(sb, gapFelToString(a));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void gapFormatPoly(StrBuffer_t* sb, const Poly_t *pol)
{
   sbAppend(sb, "[");
   for (int i = 0; i < pol->degree; ++i) {
      sbPrintf(sb,"%s,",gapFelToString(pol->data[i]));
   }
   sbPrintf(sb,"%s]",gapFelToString(pol->data[pol->degree]));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void gapFormatWord(StrBuffer_t* sb, const WgData_t *b, uint32_t n)
{
   wgDescribeWord((WgData_t *)b, n);
   int *x;
   sbAppend(sb, "[");
   for (x = b->Description; *x != -1;) {
      int first = 1;
      if (x != b->Description) {
         sbAppend(sb, ",");
      }
      sbAppend(sb, "[");
      do {
         long gen = *x++;
         if (!first) { sbAppend(sb, ",");}
         first = 0;
         sbPrintf(sb,"%ld", gen + 1);
      } while (*x != -1);
      sbAppend(sb, "]");
      ++x;
   }
   sbAppend(sb, "]");
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
