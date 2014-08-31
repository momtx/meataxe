#!/bin/sh

ME=`basename $0`
MTXBIN=`dirname $0`


if [ -z "$1" ]; then 
    echo "Usage: $ME <Command> <Args>"
    exit 1
fi

case "$1" in
    pwr*) exec $MTXBIN/zpo "$2" "$1" $3
esac

echo "zsm: unknown command $1"
exit 1
