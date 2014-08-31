#!/bin/sh
################################################################################
# MeatAxe script
################################################################################

cp null new

for i in $1.* ; do
    echo "MULTIPLYING WITH $i... "
    zmu rest $i h$i
    zcl space h$i r0 h0
    rm h$i h0

    echo "NEW FROM $i HAS "
    zef r0 r$i
    rm r0 
    zpt h0 space r$i
    mv h0 space 

    zpt r0 new r$i
    mv r0 new   
done

echo " "
echo "NEW HAS "
zef new rest

cp space space$2
cp rest rest$2
echo "RESULT SAVED TO space$2 "
echo " "

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program spin "SPIN (script)"
!synopsis 
!seealso spininit
!description
-------------------------------------------------------------------------- }}}*/
