#!/bin/sh
################################################################################
# MeatAxe script
################################################################################


set -x
ztm $1.1 $2.1 rest h1
ztm $1.2 $2.2 rest h2
ztm $1.3 $2.3 rest h3

zcl space h1 rest1 h0
zcl space h2 rest2 h0
zcl space h2 rest3 h0

zpt new rest1 rest2 rest3
zef new rest > /dev/null

zpt h0 space rest
zef h0 space

cp space space$3
cp rest rest$3

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program tcspin "TCSPIN (script)"
!synopsis 
    name
!seealso 
!description
-------------------------------------------------------------------------- }}}*/
