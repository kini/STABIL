#!/bin/sh

TESTS=*.in

make all

echo "Running tests."

for X in $TESTS ; do
    echo
    echo " ---) Testing $X ."
    echo " ---  Original STABIL:"
    ./STABIL.old "$X" | tee "${X%.in}.STABIL-old.out"
    echo " ---  STABCOL:"
    ./STABCOL "$X" | tee "${X%.in}.STABCOL.out"
    echo " ---  New STABIL:"
    ./STABIL "$X" | tee "${X%.in}.STABIL-new.out"
    echo "(---  Done testing $X ."
done
