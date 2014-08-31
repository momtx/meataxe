#!/bin/sh
################################################################################
# MeatAxe script
################################################################################


#Assumes Group,Subgroup,Vectors
#Writes conjugated permutations to z1,2 and unkondensed vectors to nsp
zpt h1 $2"1" $2"2" 
zmo h1 h2 h4
zmu h2 $1"1" h10
zmu h2 $1"2" h20
ziv h2 h3
zmu h10 h3 z1
zmu h20 h3 z2
zuk $3 h4 nsp

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program unkd "Uncondense"
!synopsis 
    name
!seealso 
!description
-------------------------------------------------------------------------- }}}*/
