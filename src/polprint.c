////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Print a polynomial
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Format a polynomial
///
/// This function converts a polynomial to a human-readable text form and appends the text to the
/// given string buffer.

void polFormat(StrBuffer* sb, const Poly_t *p)
{
   polValidate(MTX_HERE, p);
   ffSetField(p->field);
   if (p->degree == -1) {
      sbPrintf(sb, "0x^0");
   }
   const char* plus = "";
   for (int i = p->degree; i >= 0; i--) {
      if (p->data[i] != FF_ZERO) {
         sbAppend(sb, plus);
         if ((p->data[i] != FF_ONE) || (i == 0)) {
            sbPrintf(sb, "%d",ffToInt(p->data[i]));
         }
         if (i > 1) {
            sbPrintf(sb, "x^%d",i);
         } else if (i == 1) {
            sbAppend(sb, "x");
         }
         plus = "+";
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Print a polynomial
///
/// This function prints a polynomial on the standard output in a human-readable form (see
/// @ref polFormat).
/// If @a name is not 0, the name followed by an equal sign is printed before the polynomial.
/// For example, the statement <tt>polPrint("P",P)</tt> could
/// produce the following output:
/// <pre>
/// P=3x^2+x+1</pre>
///
/// @param name Name to print before the polynomial or 0.
/// @param p Pointer to the polynomial.

void polPrint(char *name, const Poly_t *p)
{
   polValidate(MTX_HERE, p);
   if (name != NULL) {
      printf("%s=",name);
   }
   StrBuffer* sb = sbAlloc(30);
   polFormat(sb, p);
   fputs(sbData(sb), stdout);
   sbFree(sb);
   if (name != NULL) {
      printf("\n");
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the human-readable form (as defined by @ref polFormat) of a polynomial in an ephemeral
/// string see @ref sbToEphemeralString).

char* polToEphemeralString(const Poly_t *p)
{
   StrBuffer* sb = sbAlloc(100);
   polFormat(sb, p);
   return sbToEphemeralString(sb);
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
