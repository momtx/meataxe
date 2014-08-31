#!/bin/sh
################################################################################
# MeatAxe script
################################################################################


#Assumes a p-element
#Writes Jordan basis to jbas and basis of endom. ring to endobas
cp $1 G1
zjo > h0
rm G1
mv P2 jbas
zmu jbas z1 h1
zmu jbas z2 h2
ziv jbas h3
zmu h1 h3 h4
zmu h2 h3 h5
zer h4 h5 < h0 > /dev/null
mv P1 endobas

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program endo "ENDO (script)"
!synopsis 
    name
!seealso 
!description
-------------------------------------------------------------------------- }}}*/
