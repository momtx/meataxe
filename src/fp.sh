#!/bin/sh
################################################################################
# MeatAxe script                                    
################################################################################

zmu z1 z2 h3 || exit 1
zad z1 z2 h4 || exit 1
zad h3 h4 h5 || exit 1
znu h5 nsp || exit 1
cp h5 eff || exit 1
#
zmu h3 z2 h6 || exit 1
zad h5 h6 h7 || exit 1
znu h7 nsp || exit 1
cp h7 eff || exit 1
#
zmu z2 h7 h8 || exit 1
zad z1 h8 h9 || exit 1
znu h9 nsp || exit 1
cp h9 eff || exit 1
#
zad z2 h9 h10 || exit 1
znu h10 nsp  || exit 1
cp h10 eff || exit 1
#
zad h3 h10 h11 || exit 1
znu h11 nsp || exit 1
cp h11 eff || exit 1
#
zad z1 h11 h12 || exit 1
znu h12 nsp || exit 1
cp h12 eff || exit 1

exit 0


/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program fp "Fingerprint"
!synopsis
    fp
!seealso @fingerprint
!description
    This program expects two matrices in `z1" and `z2', and calculates
    the fingerprint, i.e., the nullities of 6 standard words.
-------------------------------------------------------------------------- }}}*/
