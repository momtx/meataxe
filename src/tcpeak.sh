#!/bin/sh


LOG=tcpeak.log
date > $LOG

echo '===== Step 1: Chopping ====='
chop -g $6 $3 >> $LOG
chop -g $6 $4 >> $LOG

echo '===== Step 2: Find peak words ====='
pwkond -n -t $3 >> $LOG
pwkond -n -t $4 >> $LOG

echo '===== Step 3: Precondensation ====='
#mktree -g $6 $3 ${3}X${4} >> $LOG
precond ${3}X${4} $3 $4 >> $LOG

echo '===== Step 4: Transform to symmetry basis ====='
basbypeak $3 bas$3 >> $LOG
ziv bas$3 invbas$3 >> $LOG

for i in $1.[0-9] $1.[0-9][0-9] ; do
    if [ ! -r $i ] ; then continue ; fi
    zmu bas$3 $i xxx
    zmu xxx invbas$3 st$i
done

basbypeak $4 bas$4 >> $LOG
ziv bas$4 invbas$4 >> $LOG

for i in $2.[0-9] $2.[1-9][0-9] ; do
    if [ ! -r $i ] ; then continue ; fi
    zmu bas$4 $i xxx
    zmu xxx invbas$4 st$i
done

rm -f xxx

echo "===== Step 5: Condensing $1 x $2 -> result ====="
tkond -g $5 st$1 st$2 $3 $4 ${3}X${4} result >> $LOG

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program tcpeak "Tensor product condensation"
!synopsis 
    tcpeak <A> <B> <AK> <BK> <NGen> <NGenK>
!arg <A>
    Name of first representation.
!arg <B>
    Name of second representation.
!arg <AK>
    Name of first representation restricted to condensation subgroup.
!arg <BK>
    Name of second representation restricted to condensation subgroup.
!arg <NGen>
    Number of generators for representations <A> and <B>.
!arg <NGenK>
    Number of generators for condensed representations <AK> and <BK>.
!seealso tkond
!description
    This script calculates the condensed tensor product of two 
    representations. Both representations, and their condensations,
    must be available. The result is written to |result|.
-------------------------------------------------------------------------- }}}*/
