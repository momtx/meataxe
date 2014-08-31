#!/bin/sh
################################################################################
# MeatAxe script
################################################################################

ME=`basename $0`


if [ -z "$1" ]; then 
    echo "Usage: $ME <Module>"
   exit 1
fi

while  [ -z "$1" ]; do
    echo "Cleaning up $1"
    rm -f $MODULE[0-9]* $1.ssb $1.v
    shift
done

exit 0



/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program mtxclean "Clean up files"
!synopsis 
    mtxclean <Module> ...
!arg <Module>
    Name of the representation.
!description
    This program deletes all files belonging to a given module, except the
    generators.
-------------------------------------------------------------------------- }}}*/
