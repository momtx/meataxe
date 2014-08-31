#!/bin/sh
################################################################################
# MeatAxe script
################################################################################


#Writes standard basis to sbas 
zsb $1"1" $1"2" nsp sbas
ziv sbas h0
zmu sbas $1"1" h1
zmu h1 h0 z1
zmu sbas $1"2" h2
zmu h2 h0 z2

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program sb "SB (script)"
!synopsis 
    name
!seealso 
!description
-------------------------------------------------------------------------- }}}*/
