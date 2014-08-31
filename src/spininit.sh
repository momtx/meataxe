#!/bin/sh
################################################################################
# MeatAxe script
################################################################################


## Initialize "spin" with nullspace in $1 and the empty matrix in "null".

zct 0 $1 null

echo "SEED HAS "
zef $1 space

cp space rest

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program spininit "SPININIT (script)"
!synopsis 
    name
!seealso spin
!description
-------------------------------------------------------------------------- }}}*/
