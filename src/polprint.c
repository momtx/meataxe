////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Print a polynomial
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Print a polynomial
/// This function prints a polynomial on the standard output in a human-readable
/// form. If @em name is not 0, the name followed by an equal sign is printed
/// before the polynomial. For example, the statement <tt>polPrint("P",P)</tt> could
/// produce the following output:
/// <pre>
/// P=3x^2+x+1</pre>
/// @param name Name to print before the polynomial or 0.
/// @param p Pointer to the polynomial.
/// @return 0 on success, -1 on error.

void polPrint(char *name, const Poly_t *p)
{
   int i,flag = 0;

   polValidate(MTX_HERE, p);
   if (name != NULL) {
      printf("%s=",name);
   }
   ffSetField(p->field);
   if (p->degree == -1) {
      printf("0x^0");
   }
   for (i = p->degree; i >= 0; i--) {
      if (p->data[i] != FF_ZERO) {
         if (flag) {
            printf("+");
         }
         if ((p->data[i] != FF_ONE) || (i == 0)) {
            printf("%d",ffToInt(p->data[i]));
         }
         if (i > 1) {
            printf("x^%d",i);
         } else if (i == 1) {
            printf("x");
         }
         flag = 1;
      }
   }
   if (name != NULL) {
      printf("\n");
   }
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
