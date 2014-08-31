#!/bin/sh
################################################################################
# MeatAxe script
################################################################################

#!/bin/sh

## Split by hand

echo "SUBSPACE HAS "
zef $2 spc

echo "NEXT RANKS SHOULD BE 0... "
for i in $1.* ; do
    zmu spc $i h$i
    zcl spc h$i h0 s$i
    
    echo "FOR $i: "
    zef h0 hh0
    rm h$i h0 hh0

    zqt -i spc $i q$i
done


exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program sp "SP (script)"
!synopsis 
    name
!seealso 
!description
-------------------------------------------------------------------------- }}}*/
