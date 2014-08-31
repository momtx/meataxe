#include "meataxe.h"
#include "config.h"

#define MKSTRING2(x) #x
#define MKSTRING(x) MKSTRING2(x)

static char Version_[] = 
   "@(#)$MeatAxeVersion: "
   MTX_VERSION "." MTX_BUILD
   " ZZZ=" MKSTRING(ZZZ)
   MTX_CONFIG_STRING
   " " __DATE__ " " __TIME__
   " $";
char *MtxVersion = Version_ + 21;
