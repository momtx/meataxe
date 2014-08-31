#!/bin/sh
################################################################################
# MeatAxe script
################################################################################


#Assumes word,group1,group2
#Writes standardized gen. to s11,2 and s21,2
zsm mw$1 $2"1" $2"2" h1 h2
zsb $2"1" $2"2" h2 bas1
ziv bas1 h0
zmu bas1 $2"1" h1
zmu h1 h0 s11
zmu bas1 $2"2" h2
zmu h2 h0 s12
zsm mw$1 $3"1" $3"2" h1 h2
zsb $3"1" $3"2" h2 bas2
ziv bas2 h0
zmu bas2 $3"1" h1
zmu h1 h0 s21
zmu bas2 $3"2" h2
zmu h2 h0 s22
cmp s11 s21
cmp s12 s22

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program eq "EQ (script)"
!synopsis 
    name
!seealso 
!description
-------------------------------------------------------------------------- }}}*/
