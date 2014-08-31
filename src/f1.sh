#!/bin/sh
################################################################################
# MeatAxe script
################################################################################

zmu z1 z2 h3 || exit 1
zad z1 z2 h4 || exit 1
zad h3 h4 h5 || exit 1
znu h5 nsp   || exit 1
cp h5 eff
exit 0
 

/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program f1 "First fingerprint"
!synopsis 
    f1
!seealso @fingerprint
!description
    This program expects two matrices in `z1' and `z2', calculates
    the first standard word (A+B+AB) in `eff' and its null-space
    in `nsp'.
-------------------------------------------------------------------------- }}}*/
