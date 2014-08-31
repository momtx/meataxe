#!/bin/sh
################################################################################
# MeatAxe script - see end of file for documentation
################################################################################

ME=`basename $0`
MODULE=$1
if [ ! -r $1.cfinfo ] ; then
    echo "$ME: $1.cfinfo not found -- run 'chop $1'"
    exit 1
fi

echo "Calculating peak words for $MODULE"
pwkond -Qt $MODULE || exit 1
echo "Calculating the endomorphism ring $MODULE.endo and $MODULE.endo.lrr"
mkhom -Q -r 1 $MODULE $MODULE $MODULE.endo || exit 1
echo "Calculating peak words of $MODULE.endo.lrr"
chop -Q -i $MODULE.endo.lrr || exit 1
pwkond -Qt $MODULE.endo.lrr || exit 1
echo "Calculating the socle of $MODULE.endo.lrr"
soc -l 1 -Q $MODULE.endo.lrr || exit 1
echo "Decomposing $MODULE"
decomp $MODULE $MODULE.endo || exit 1
exit 0



/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program dec "Decompose a module"
!synopsis 
    dec <Module>
!arg <Module>
    Name of the representation.
!arg <End>
    Name for the endomorphism ring.
!seealso chop decomp mkhom soc pwkond
!index "decomposing a module"
!description
    This program finds the indecomposable direct summands of a module.
    `<Module>' is the name of the module. The program calculates the
    endomporphism ring of the module as `<Module>.endo'. It assumes
    than the module has been chopped.
-------------------------------------------------------------------------- }}}*/
